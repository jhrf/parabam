import pysam
import parabam
import shutil
import time
import sys
import os

from multiprocessing import Queue

from parabam.interface import util as ut
from parabam.interface import merger


class TaskSubset(parabam.Task):
	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(TaskSubset, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
		self._merge_types = const.merge_types
		
	def produceResultsDict(self):
		master = pysam.Samfile(self._master_file_path,"rb")

		temp_file_names = {}
		temp_file_objects = {}
		counts = {}

		for mT in self._merge_types:
			temp_file_names[mT] = self.__mkTmpName__(mT)
			temp_file_objects[mT]  = pysam.Samfile(temp_file_names[mT],"wb",template=master)
			counts[mT] = 0

		engine = self.const.user_engine
		user_constants = self.const.user_constants

		for alig in iter(self._taskSet):
			#If true, we will add this read to the output file
			if engine(alig,user_constants,master):
				#As this is a simpflified subsetter, we only have one merge_type
				temp_file_objects[self._merge_types[0]].write(alig) #Write telomeric read to file
				counts[self._merge_types[0]] += 1

		master.close() #Close the master bamfile
		map(lambda (mT,fil): fil.close(),temp_file_objects.items()) #close all the other bams

		results = {}
		results["source"] = self._source
		results["filnms"] = temp_file_names
		results["counts"] = counts

		return results
		
class HandlerSubset(parabam.Handler):

	def __init__(self,inqu,outqu,const):
		super(HandlerSubset,self).__init__(inqu,const,maxDestroy=len(const.sources))

		self._sources = const.sources
		self._merge_types = const.merge_types

		#Setup stores and connection to merge proc
		self._mrgStores = {}
		for src in self._sources:
			self._mrgStores[src] = {} 
			for mT in self._merge_types:
				self._mrgStores[src][mT] = []

		self._mergequeue = outqu
		self._mergecount = 0

	def newResultAction(self,newResult,**kwargs):
		resDict = newResult.results
		source = resDict["source"]
		self.__autoHandle__(resDict,source)

		for mT in self._merge_types:
			if resDict["counts"][mT] > 0:
				self._mrgStores[source][mT].append(resDict["filnms"][mT])
			else:
				os.remove(resDict["filnms"][mT])

	def periodicAction(self,iterations):
		for src in self._sources:
			for mT in self._merge_types:
				self.__testMergeStore__(self._filenameOut[src][mT],self._mrgStores[src][mT],src,mT)

	def __testMergeStore__(self,outnm,store,src,mT):
		if len(store) > 30:
			sys.stdout.flush()
			self.__addMergeTask__(name=outnm,results=store,merge_type=mT,source=src,total=self._stats[src]["total"])
			self._mergecount += 1
			#Remove the temp file which has been merged
			store[:] = [] #Wipe the store clean, these have been merged
			sys.stdout.flush()

	def __addMergeTask__(self,name,results,merge_type,source,total,destroy=False):
		res = merger.ResultsMerge(name=name,results=list(results),
							merge_type=merge_type,source=source,destroy=destroy,total=total)
		self._mergequeue.put(res)

	def __totalSum__(self):
		return sum(map(lambda s : self._stats[s]["total"],self._sources))

	def handlerExitFunc(self,**kwargs):
		self._updateFunc("[Result] Processed %d reads from bam files\n" % (self.__totalSum__(),))
		self._updateFunc("[Status] Waiting for merge operation to finish...\n")

		last_source = len(self._sources) - 1
		last_mrg = len(self._merge_types) - 1

		for src_count,src in enumerate(self._sources):
			for merge_count,mT in enumerate(self._merge_types):
				destroy = True if (src_count == last_source and merge_count == last_mrg) else False
				self.__addMergeTask__(name=self._filenameOut[src][mT],
									results=self._mrgStores[src][mT],mrgTyp=mT,
									source=src,total=self._stats[src]["total"],
									destroy=destroy)

class ProcessorSubset(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args):
		super(ProcessorSubset,self).__init__(outqu,const,TaskClass,task_args)
		self._source = task_args[0] #Defined in the run function within Interface

	def __getMasterBam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def addToCollection(self,master,alig,collection):
		collection.append(alig)

	def preProcActivity(self,master_file_path):
		pass

class Interface(parabam.Interface):

	def __init__(self,tempDir,exe_dir):
		super(Interface,self).__init__(tempDir,exe_dir)
	
	def run_cmd(self,parser):

		cmd_args = parser.parse_args()

		module = __import__(cmd_args.instruc, fromlist=[''])
		user_engine = module.engine
		user_constants = {}
		module.set_constants(user_constants)

		self.run(
			input_bams=cmd_args.input,
			outputs= cmd_args.output,
			proc= cmd_args.p,
			chunk= cmd_args.c,
			verbose= cmd_args.v,
			merge_types= ["subset"],
			user_constants = user_constants,
			user_engine = user_engine
			)
	
	def run(self,input_bams,outputs,proc,chunk,verbose,merge_types,user_constants,user_engine):

		if outputs == None or len(outputs) == len(input_bams): 
			for input_group,output_group in self.__getGroup__(input_bams,outputs):
				
				quPrim = Queue()
				quMerg = Queue()

				outFiles = dict([(src,{}) for src in output_group])
				master_file_path = {}

				for mst,src in zip(input_group,output_group):
					master_file_path[src] = mst
					if len(merge_types) > 1: 
						for typ in merge_types:
							outFiles[src][typ] = "%s/%s_%s.bam" % (self._tempDir,src.replace(".bam",""),typ,)
					else:
						outFiles[src][merge_types[0]] = "%s/%s" % (self._tempDir,src,)
						#This has been removed and the merger now looks after output
						#ut.create_subset_bam(mst,outFiles[src][merge_types[0]])
						
				if verbose: self.__reportFileNames__(outFiles)

				procrs = []
				handls = []

				const = parabam.Const(outFiles=outFiles,
									tempDir=self._tempDir,
									master_file_path=master_file_path,
									chunk=chunk,proc=(proc // len(input_group)),
									verbose=verbose,thresh=0,
									merge_types=merge_types,
									sources=output_group,
									exe_dir=self._exe_dir,
									user_constants=user_constants,
									user_engine=user_engine)

				for src in output_group:
					procrs.append(ProcessorSubset(outqu=quPrim,
												const=const,
												TaskClass=TaskSubset,
												task_args=[src]))

				handls.append(HandlerSubset(inqu=quPrim,outqu=quMerg,const=const))
				handls.append(merger.HandlerMerge(inqu=quMerg,const=const))

				lev = parabam.Leviathon(procrs,handls,100000)
				lev.run()

				#Move the complete telbams out of the tempdir to the working dir
				#Only do this if we custom generated the file locations.
				self.__moveOutputFiles__(outFiles)
		else:
			print "[Error] Please provide the same amount of output names as input files"
			print "[Error] `parabam subset` cannot continue"

	def __reportFileNames__(self,outFiles):
		print "[Update] parabam will output the following files:"
		for src,mrgTypDict in outFiles.items():
			for mrgTyp,outFilePath in mrgTypDict.items():
				print "\t%s" % (outFilePath.split("/")[-1],)
		print ""

	def __moveOutputFiles__(self,outFiles):
		for src, mrgTypDict in outFiles.items():
			for mrgTyp,outFilePath in mrgTypDict.items():
				try:
					shutil.move(outFilePath,"./") #./ being the current working dir
				except shutil.Error,e:
					alt_filnm = "./%s_%s_%d.bam" % (src,mrgTyp,time.time()) 
					print "[Error] Output file may already exist, you may not" \
					"have correct permissions for this file"
					print "[Status]Trying to create using unique filename:"
					print "\t\t%s" % (alt_filnm,)
					shutil.move(outFilePath,alt_filnm)

	def __getGroup__(self,bams,names,multi=2):
		for i in xrange(0,len(bams),multi):
			yield (bams[i:i+multi],names[i:i+multi])

	def getParser(self):
		#argparse imported in ./interface/parabam 
		parser = ut.default_parser()

		# parser.add_argument('-t',type=int,metavar="INT",nargs='?',default=2
		# 	,help="How many TTAGGG sequences you wish to observe before" \
		# 	" accepting a sequence as telomeric.")

		return parser 