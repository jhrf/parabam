import pysam
import parabam
import shutil
import time
import sys
import os

from multiprocessing import Queue

from parabam.interface import util as ut
from parabam.interface import merger

from abc import ABCMeta, abstractmethod

class TaskStat(parabam.Task):

	__metaclass__ = ABCMeta

	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(TaskStat, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
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

		master.close() #Close the master bamfile
		map(lambda (subset,fil): fil.close(),temp_file_objects.items()) #close all the other bams

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
		
class HandlerStat(parabam.Handler):

	def __init__(self,inqu,const):
		super(HandlerStat,self).__init__(inqu,const,maxDestroy=len(const.sources))

		self._sources = const.sources
		self._mergecount = 0

	def newResultAction(self,newResult,**kwargs):
		resDict = newResult.results
		source = resDict["source"]
		self.__autoHandle__(resDict,source)



	def periodicAction(self,iterations):
		pass

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

class ProcessorStat(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args):
		super(ProcessorStat,self).__init__(outqu,const,TaskClass,task_args)
		self._source = task_args[0] #Defined in the run function within Interface

	def __getMasterBam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def addToCollection(self,master,alig,collection):
		collection.append(alig)

	def preProcActivity(self,master_file_path):
		pass

class UserStructure(object):
	def __init__(self,struc_type,record_type,data):
		self.struc_type = struc_type
		self.record_type = record_type 
		self.data = data
		self.enter_data = None

	def max_decision(self,new,existing):
		return new > existing

	def min_decision(self,new,exisiting):
		return new < existing

class NumericStructure(UserStructure):
	def __init__(self,struc_type,record_type,data,log_scaling=False):
		super(NumericStructure,self).__init__(struc_type,record_type,data)
		self.log_scaling = log_scaling
		if record_type == "max":
			self.enter_data = numeric_add_max
		elif record_type == "max":
			self.enter_data = numeric_add_min
		else:
			self.enter_data = numeric_add_cumu

	def numeric_add_cumu(self,new):
		data += new

	def numeric_add_max(self,new):
		if self.max_decision(new,data):
			data = new

	def numeric_add_min(self,new):
		if self.min_decision(new,data):
			data = new

class ArrayStructure(UserStructure):
	def __init__(self,struc_type,record_type,data):
		super(NumericStructure,self).__init__(struc_type,record_type,data)
		self.log_scaling = log_scaling

		if record_type == "max":
			self.enter_data = numeric_add_max
		elif record_type == "max":
			self.enter_data = numeric_add_min
		else:
			self.enter_data = numeric_add_cumu

	def array_add_max(self,new,coords):
		existing = data[coords]
		if self.max_decision(new,existing):
			data[coords] = new

	def array_add_min(self,new,coords):
		existing = data[coords]
		if self.min_decision(new,existing):
			data[coords] = new

	def array_add_cumu(self,new,coords):
		data[coords] += new


class Interface(parabam.Interface):

	def __init__(self,tempDir,exe_dir):
		super(Interface,self).__init__(tempDir,exe_dir)
	
	def run_cmd(self,parser):

		cmd_args = parser.parse_args()

		verbose = cmd_args.v
		module = __import__(cmd_args.instruc, fromlist=[''])
		user_engine = module.engine
		user_constants = {}
		user_structures = {}
		module.set_constants(user_constants)
		module.set_structures(user_structures)

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
			user_constants = user_constants,
			user_engine = user_engine,
			user_structures = user_structures
			engine_is_class = False
			)
	
	def run(self,input_bams,outputs,proc,chunk,verbose,user_constants,
			user_engine,engine_is_class,user_structures,multi=4):

		if not outputs or not len(outputs) == len(input_bams):
			print "[Warning] Output files will use automatic default naming scheme \n"\
			"\t\tTo specify output names ensure an output name is provided for each input BAM"
			outputs = [ ut.get_bam_basename(b) for b in input_bams ]

		#AT SOME POINT WE SHOULD HANDLE UNSORTED BAMS. EITHER HERE OR AT THE PROCESSOR

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
								sources=output_group,
								exe_dir=self._exe_dir,
								user_constants=user_constants,
								user_structures=user_structures,
								user_engine=user_engine)

			for src in output_group:
				if engine_is_class:
					if not issubclass(user_engine,TaskStat):
						raise Exception("[ERROR]\tThe class provided to parabam multiset must be a subclass of\n"\
										"\tparabam.interface.subset.TaskStat. Please consult the parabam manual.")
					cur_args = [src]
					procrs.append(ProcessorStat(outqu=quPrim,
											const=const,
											TaskClass=user_engine,
											task_args=cur_args))
				else:
					if len(subset_types) == 1:
						run_class = SingleSet
					else:
						run_class = MultiSet
					procrs.append(ProcessorStat(outqu=quPrim,
												const=const,
												TaskClass=run_class,
												task_args=[src]))

			handls.append(HandlerStat(inqu=quPrim,const=const))

			lev = parabam.Leviathon(procrs,handls,100000)
			lev.run()

			#Move the complete telbams out of the tempdir to the working dir
			#Only do this if we custom generated the file locations.
			self.__moveOutputFiles__(outFiles)

	def getParser(self):
		#argparse imported in ./interface/parabam 
		parser = ut.default_parser()

		# parser.add_argument('-t',type=int,metavar="INT",nargs='?',default=2
		# 	,help="How many TTAGGG sequences you wish to observe before" \
		# 	" accepting a sequence as telomeric.")

		return parser 