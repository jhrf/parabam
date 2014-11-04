import pysam
import parabam
import time
import sys
import os
import gc

from multiprocessing import Queue

from parabam.interface import util as ut
from parabam.interface import merger

from abc import ABCMeta, abstractmethod

class TaskSubset(parabam.Task):

	__metaclass__ = ABCMeta

	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(TaskSubset, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
		self._subset_types = const.subset_types
		self._counts = {}
		self._temp_file_names = {}
		self._temp_file_objects = {}

	def produceResultsDict(self):
		master = pysam.Samfile(self._master_file_path,"rb")

		temp_file_names = self._temp_file_names
		temp_file_objects = self._temp_file_objects
		counts = self._counts

		for subset in self._subset_types:
			temp_file_names[subset] = self.__mkTmpName__(subset)
			temp_file_objects[subset]  = pysam.Samfile(temp_file_names[subset],"wb",template=master)
			counts[subset] = 0

		self.handle_task_set(self._taskSet,master)

		#Close the master bamfile
		for subset,temp_obj in temp_file_objects.items():
			temp_obj.close() #close all the other bams
		master.close()

		results = {}
		results["source"] = self._source
		results["filnms"] = temp_file_names
		results["counts"] = counts

		return results

	def write_to_subset_bam(self,subset_type,read):
		self._counts[subset_type] += 1
		self._temp_file_objects[subset_type].write(read)

	@abstractmethod
	def handle_task_set(self,task_set,master):
		pass

class MultiSet(TaskSubset):
	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(MultiSet, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const,
										source=source)
	
	def handle_task_set(self,task_set,master):

		engine = self.const.user_engine
		user_constants = self.const.user_constants
		subset_write = self.write_to_subset_bam
		subset_types = self.const.subset_types

		for read in task_set:
			subset_decision = engine(read,user_constants,master)
			if type(subset_decision) == int:
				if not subset_decision == -1:
					subset_write(subset_types[subset_decision],read)
			elif type(subset_decision) == tuple:
				for subset in subset_decision:
					subset_write(subset_types[subset],read)					

class SingleSet(TaskSubset):
	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(SingleSet, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const,
										source=source)
	
	def handle_task_set(self,task_set,master):
		engine = self.const.user_engine
		user_constants = self.const.user_constants
		subset_write = self.write_to_subset_bam
		subset_types = self.const.subset_types

		for read in task_set:
			if engine(read,user_constants,master):
				subset_write(subset_types[0],read)
		
class HandlerSubset(parabam.Handler):

	def __init__(self,inqu,outqu,const):
		super(HandlerSubset,self).__init__(inqu,const,maxDestroy=len(const.sources))

		self._sources = const.sources
		self._subset_types = const.subset_types

		#Setup stores and connection to merge proc
		self._mrgStores = {}
		for src in self._sources:
			self._mrgStores[src] = {} 
			for subset in self._subset_types:
				self._mrgStores[src][subset] = []

		self._mergequeue = outqu
		self._mergecount = 0

	def newResultAction(self,newResult,**kwargs):
		resDict = newResult.results
		source = resDict["source"]
		self.__autoHandle__(resDict,source)

		for subset in self._subset_types:
			if resDict["counts"][subset] > 0:
				self._mrgStores[source][subset].append(resDict["filnms"][subset])
			else:
				os.remove(resDict["filnms"][subset])

	def periodicAction(self,iterations):
		for src in self._sources:
			for subset in self._subset_types:
				self.__testMergeStore__(self._filenameOut[src][subset],self._mrgStores[src][subset],src,subset)

	def __testMergeStore__(self,outnm,store,src,subset):
		if len(store) > 10:
			sys.stdout.flush()
			self.__addMergeTask__(name=outnm,results=store,subset_type=subset,source=src,total=self._stats[src]["total"])
			self._mergecount += 1
			#Remove the temp file which has been merged
			store[:] = [] #Wipe the store clean, these have been merged
			sys.stdout.flush()

	def __addMergeTask__(self,name,results,subset_type,source,total,destroy=False):
		res = merger.ResultsMerge(name=name,results=list(results),
							subset_type=subset_type,source=source,destroy=destroy,total=total,time_added=time.time())
		self._mergequeue.put(res)

	def __totalSum__(self):
		return sum(map(lambda s : self._stats[s]["total"],self._sources))

	def handlerExitFunc(self,**kwargs):
		self._updateFunc("[Result] Processed %d reads from bam files\n" % (self.__totalSum__(),))
		self._updateFunc("[Status] Waiting for merge operation to finish...\n")

		last_source = len(self._sources) - 1
		last_mrg = len(self._subset_types) - 1

		for src_count,src in enumerate(self._sources):
			for merge_count,mT in enumerate(self._subset_types):
				destroy = True if (src_count == last_source and merge_count == last_mrg) else False
				self.__addMergeTask__(name=self._filenameOut[src][mT],
									results=self._mrgStores[src][mT],subset_type=mT,
									source=src,total=self._stats[src]["total"],
									destroy=destroy)

class ProcessorSubset(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args,debug=False):
		# if const.fetch_region:
		# 	debug = False 
		super(ProcessorSubset,self).__init__(outqu,const,TaskClass,task_args,debug)
		self._source = task_args[0] #Defined in the run function within Interface

	def __getMasterBam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def addToCollection(self,master,alig,collection):
		collection.append(alig)

	def __getNextAlig__(self,masterBam):
		if not self.const.fetch_region:
			for alig in masterBam.fetch(until_eof=True):
				yield alig
		else:
			for alig in masterBam.fetch(region=self.const.fetch_region):
				yield alig

	def preProcActivity(self,master_file_path):
		pass

class Interface(parabam.Interface):

	def __init__(self,tempDir,exe_dir):
		super(Interface,self).__init__(tempDir,exe_dir)
	
	def run_cmd(self,parser):

		cmd_args = parser.parse_args()

		verbose = cmd_args.v
		module = __import__(cmd_args.instruc, fromlist=[''])
		user_engine = module.engine
		user_constants = {}
		module.set_constants(user_constants)

		if hasattr(module,"get_subset_types"):
			if verbose: 
				print "[Status] Multiple subset type identified"
			subset_types = module.get_subset_types()
		else:
			subset_types = ["subset"]

		self.run(
			input_bams=cmd_args.input,
			outputs= cmd_args.output,
			proc= cmd_args.p,
			chunk= cmd_args.c,
			verbose= verbose,
			subset_types= subset_types,
			user_constants = user_constants,
			user_engine = user_engine,
			engine_is_class = False,
			fetch_region = cmd_args.region
			)
	
	def run(self,input_bams,outputs,proc,chunk,subset_types,
			user_constants,user_engine,fetch_region=None,
			multi=4,keep_in_temp=False,engine_is_class=False,verbose=False):

		if not outputs or not len(outputs) == len(input_bams):
			print "[Warning] Output files will use default naming scheme \n"\
			"\t\tTo specify output names ensure a name is provided for each input BAM"
			outputs = [ ut.get_bam_basename(b) for b in input_bams ]

		#AT SOME POINT WE SHOULD HANDLE UNSORTED BAMS. EITHER HERE OR AT THE PROCESSOR
		final_files = []

		for input_group,output_group in self.__getGroup__(input_bams,outputs,multi=multi):
			
			quPrim = Queue()
			quMerg = Queue()

			outFiles = dict([(src,{}) for src in output_group])
			master_file_path = {}

			for mst,src in zip(input_group,output_group):
				master_file_path[src] = mst
				for typ in subset_types:
					outFiles[src][typ] = "%s/%s_%s.bam" % (self._tempDir,src.replace(".bam",""),typ,)
					
			if verbose: self.__reportFileNames__(outFiles)

			procrs = []
			handls = []

			const = parabam.Const(outFiles=outFiles,
								tempDir=self._tempDir,
								master_file_path=master_file_path,
								chunk=chunk,proc=(proc // len(input_group)),
								verbose=verbose,thresh=0,
								subset_types=subset_types,
								sources=output_group,
								exe_dir=self._exe_dir,
								user_constants=user_constants,
								user_engine=user_engine,
								fetch_region=fetch_region)

			for src in output_group:
				if engine_is_class:
					if not issubclass(user_engine,TaskSubset):
						raise Exception("[ERROR]\tThe class provided to parabam multiset must be a subclass of\n"\
										"\tparabam.interface.subset.TaskSubset. Please consult the parabam manual.")
					cur_args = [src]
					procrs.append(ProcessorSubset(outqu=quPrim,
											const=const,
											TaskClass=user_engine,
											task_args=cur_args))
				else:
					if len(subset_types) == 1:
						run_class = SingleSet
					else:
						run_class = MultiSet
					procrs.append(ProcessorSubset(outqu=quPrim,
												const=const,
												TaskClass=run_class,
												task_args=[src]))

			handls.append(HandlerSubset(inqu=quPrim,outqu=quMerg,const=const))
			handls.append(merger.HandlerMerge(inqu=quMerg,const=const))

			lev = parabam.Leviathon(procrs,handls,100000)
			lev.run()
			del lev

			#Move the complete telbams out of the tempdir to the working dir
			#Only do this if we custom generated the file locations.
			if keep_in_temp:
				for source,subset_paths in outFiles.items():
					for subset,path in subset_paths.items():
						final_files.append(path)
			else:
				final_files.extend(self.__moveOutputFiles__(outFiles))
			
			gc.collect()

		return final_files

	def getParser(self):
		#argparse imported in ./interface/parabam 
		parser = ut.default_parser()

		parser.add_argument('-r','--region',type=str,metavar="REGION",nargs='?',default=None
			,help="The subset process will be run only on reads from this region. \
			Regions should be colon seperated as specifiec by samtools (eg \'chr1:1000,5000\')")

		return parser 