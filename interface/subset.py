import pysam
import parabam
import time
import sys
import os
import gc
import shutil

from multiprocessing import Queue

from parabam.interface import merger

from abc import ABCMeta, abstractmethod

class TaskSubset(parabam.Task):

	def __init__(self,task_set,outqu,curproc,destroy,const,source):
		super(TaskSubset, self).__init__(task_set=task_set,
										outqu=outqu,
										curproc=curproc*len(const.subset_types),
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
		self._subset_types = const.subset_types
		self._counts = {}
		self._temp_paths = {}
		self._temp_objects = {}

	def __generate_results__(self):
		master = pysam.Samfile(self._master_file_path,"rb")

		temp_paths = self._temp_paths
		temp_objects = self._temp_objects
		counts = self._counts

		for subset in self._subset_types:
			temp_paths[subset] = self.__get_temp_path__(subset)
			temp_objects[subset]  = pysam.Samfile(temp_paths[subset],"wb",template=master)
			counts[subset] = 0

		self.__handle_task_set__(self._task_set,master)

		#Close the master bamfile
		for subset,file_object in temp_objects.items():
			file_object.close() #close all the other bams
		del self._temp_objects
		master.close()

		results = {}
		results["source"] = self._source
		results["temp_paths"] = temp_paths
		results["counts"] = counts

		return results

	def __write_to_subset_bam__(self,subset_type,read):
		self._counts[subset_type] += 1
		self._temp_objects[subset_type].write(read)

	def __write_to_subset_bam_pair__(self,subset_type,reads):
		for read in reads:
			self._counts[subset_type] += 1
			self._temp_objects[subset_type].write(read)

	def __handle_task_set__(self,task_set,master):

		engine = self.const.user_engine
		user_constants = self.const.user_constants
		subset_types = self.const.subset_types

		subset_write = self.__write_to_subset_bam__
		if self.const.pair_process:
			subset_write = self.__write_to_subset_bam_pair__

		for read in task_set:
			subset_decision = engine(read,user_constants,master)
			
			if type(subset_decision) == bool:
				if subset_decision:
					subset_write(subset_types[0],read)			
			elif type(subset_decision) == list:
				for subset,cur_read in subset_decision:
					subset_write(subset,cur_read)
			elif type(subset_decision) == tuple:
				for subset in subset_decision:
					subset_write(subset_types[subset],read)
			elif type(subset_decision) == int:
				if not subset_decision == -1:
					subset_write(subset_types[subset_decision],read)
		
class HandlerSubset(parabam.Handler):

	def __init__(self,inqu,outqu,const,destroy_limit):
		super(HandlerSubset,self).__init__(inqu,const,destroy_limit=destroy_limit)

		self._sources = const.sources
		self._subset_types = const.subset_types

		#Setup stores and connection to merge proc
		self._merge_stores = {}
		for source in self._sources:
			self._merge_stores[source] = {} 
			for subset in self._subset_types:
				self._merge_stores[source][subset] = []

		self._mergequeue = outqu
		self._mergecount = 0

	def __get_total_processed_reads__(self):
		total = 0
		for source,deep_stat in self._stats.items():
			for name,value in deep_stat.items():
				if name == "total":
					total += value
		return total

	def __new_package_action__(self,new_package,**kwargs):
		results = new_package.results
		source = results["source"]
		self.__auto_handle__(results,source)

		for subset in self._subset_types:
			self._merge_stores[source][subset].append((results["counts"][subset],results["temp_paths"][subset],))

	def __periodic_action__(self,iterations):
		for source in self._sources:
			for subset in self._subset_types:
				self.__test_merge_store__(self._output_paths[source][subset],self._merge_stores[source][subset],source,subset)

	def __test_merge_store__(self,outnm,store,source,subset):
		if len(store) > 1:
			self.__add_merge_task__(name=outnm,results=store,subset_type=subset,source=source,total=self._stats[source]["total"])
			self._mergecount += 1
			#Remove the temp file which has been merged
			store[:] = [] #Wipe the store clean, these have been merged
			gc.collect()

	def __add_merge_task__(self,name,results,subset_type,source,total,destroy=False):
		res = merger.MergePackage(name=name,results=list(results),
							subset_type=subset_type,source=source,destroy=destroy,total=total,time_added=time.time())
		self._mergequeue.put(res)

	def __total_reads__(self):
		return sum(map(lambda s : self._stats[s]["total"],self._sources))

	def __handler_exit__(self,**kwargs):
		if self._verbose:
			self.__standard_output__("\n[Status] Processed %d reads from bam files\n" % (self.__total_reads__(),))
			self.__standard_output__("[Status] Waiting for merge operation to finish...\n")

		for source in self._sources:
			for subset in self._subset_types:
				self.__add_merge_task__(name=self._output_paths[source][subset],
								results=self._merge_stores[source][subset],subset_type=subset,
								source=source,total=self._stats[source]["total"],
								destroy=True)

class ProcessorSubset(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args,debug=False):
		# if const.fetch_region:
		# 	debug = False 
		super(ProcessorSubset,self).__init__(outqu,const,TaskClass,task_args,debug)
		self._source = task_args[0] #Defined in the run function within Interface

	def __get_master_bam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def __add_to_collection__(self,master,alig,collection):
		collection.append(alig)

	def __get_next_alig__(self,master_bam):
		if not self.const.fetch_region:
			for alig in master_bam.fetch(until_eof=True):
				yield alig
		else:
			for alig in master_bam.fetch(region=self.const.fetch_region):
				yield alig

	def __pre_processor__(self,master_file_path):
		pass

class PairProcessor(ProcessorSubset):
	def __init__(self,outqu,const,TaskClass,task_args,debug=False):
		super(PairProcessor,self).__init__(outqu,const,TaskClass,task_args,debug=False)
		self._loners = {}
		self._loner_count = 0
		master = pysam.Samfile(self._master_file_path[task_args[0]],"rb")
		self._loners_object = pysam.Samfile("%s/loners" % (self._temp_dir,),"wb",template=master)
		master.close()

	def __add_to_collection__(self,master,item,collection):
		loners = self._loners

		try:
			mate = loners[item.qname] 
			del loners[item.qname]
			self._loner_count -= 1
			collection.append( (item,mate,) )
		except KeyError:
			#Could implement a system where by long standing
			#unpaired reads are stored to be run at the end 
			#of the program, otherwise we risk clogging memory
			loners[item.qname] = item
			self._loner_count += 1

		sys.stdout.write("\r %d" % (self._loner_count,))

		if self._loner_count > 200000:
		 	self._loner_count = 0
		 	self.__stash_loners__(loners)
		 	del self._loners
		 	gc.collect()
		 	self._loners = {}

	def __stash_loners__(self,loners):
		for name,read in loners.items():
			self._loners_object.write(read)

	def __end_processing__(self,master):
		super(PairProcessor,self).__end_processing__(master)
		self.__stash_loners__(self._loners)
		self._loners_object.close()

class Interface(parabam.UserInterface):
	"""The interface to parabam subset.
	Users will primarily make use of the ``run`` function."""

	def __init__(self,temp_dir,exe_dir):
		super(Interface,self).__init__(temp_dir,exe_dir)

	def run_cmd(self,parser):
		cmd_args = parser.parse_args()

		verbose = cmd_args.v

		module,user_engine,user_constants = self.__get_module_and_vitals__(cmd_args.instruc)

		if hasattr(module,"get_subset_types"):
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
			fetch_region = cmd_args.region,
			pair_process=cmd_args.m,
			side_by_side = cmd_args.s,
			debug = cmd_args.debug
			)
	
	def run(self,input_bams,outputs,proc,chunk,subset_types,
			user_constants,user_engine,fetch_region=None,side_by_side=2,
			keep_in_temp=False,engine_is_class=False,verbose=False,
			pair_process=False,debug=False):

		"""This is function is invoked in order to run parabam subset programatically.

		Arguments:
			input_bams (list): BAM Files we wish to analyse
			outputs (list): Names for the output files. If not specified <bam_name>_<subset_types>.bam convention will be used
			proc (int): The maximum number of tasks that will be run at one time
			chunk (int): The amount of reads that a task will hold
			subset_types (list): A list of names which correspond to the subset types
			user_constants (dict): A dictionary initialised with object referenced by the user engine
			user_engine (function): A function or class to be run on each read
			verbose (int): Expects an int from 0 to 2. The level of output produced by parabam
			keep_in_temp (bool): Files will be kept in temp file after processing. Useful for incorporation into pipelines
			side_by_side (int): The amount of BAMs to be run in side-by-side mode. Further concurrency
			pair_process (bool): Pair processing mode will be initialised (feature in beta)
			debug (bool): subset will only process the first 5million reads from the BAM files."""

		if not outputs or not len(outputs) == len(input_bams):
			print "[Status] Using default naming scheme."
			outputs = [ self.__get_basename__(b) for b in input_bams ]

		#AT SOME POINT WE SHOULD HANDLE UNSORTED BAMS. EITHER HERE OR AT THE PROCESSOR
		final_files = []

		for input_group,output_group in self.__get_group__(input_bams,outputs,multi=side_by_side):
			
			task_qu = Queue()
			merge_qu = Queue()

			output_paths = dict([(source,{}) for source in output_group])
			master_file_path = {}

			for mst,source in zip(input_group,output_group):
				master_file_path[source] = mst
				for typ in subset_types:
					output_paths[source][typ] = "%s/%s_%s.bam" % (self._temp_dir,source.replace(".bam",""),typ,)
					
			if verbose: self.__report_file_names__(output_paths)

			procrs = []
			handls = []

			const = parabam.Const(output_paths=output_paths,
								temp_dir=self._temp_dir,
								master_file_path=master_file_path,
								chunk=chunk,proc=(proc // len(input_group)),
								verbose=verbose,thresh=0,
								subset_types=subset_types,
								sources=output_group,
								exe_dir=self._exe_dir,
								user_constants=user_constants,
								user_engine=user_engine,
								fetch_region=fetch_region,
								pair_process=pair_process)

			for source in output_group:
				processor_class = ProcessorSubset
				if pair_process:
					processor_class = PairProcessor

				if engine_is_class:
					if not issubclass(user_engine,TaskSubset):
						raise Exception("[ERROR]\tThe class provided to parabam multiset must be a subclass of\n"\
										"\tparabam.interface.subset.TaskSubset. Please consult the parabam manual.")
					cur_args = [source]
					procrs.append(processor_class(outqu=task_qu,
											const=const,
											TaskClass=user_engine,
											task_args=cur_args,
											debug=debug))
				else:
					procrs.append(processor_class(outqu=task_qu,
												const=const,
												TaskClass=TaskSubset,
												task_args=[source],
												debug = debug))

			destroy_limit=len(const.sources)
			handls.append(HandlerSubset(inqu=task_qu,outqu=merge_qu,const=const,destroy_limit=destroy_limit))
			handls.append(merger.HandlerMerge(inqu=merge_qu,const=const,destroy_limit=destroy_limit))

			if verbose == 1: 
				update_interval = 100
			else:
				update_interval = 1

			lev = parabam.Leviathon(procrs,handls,update_interval)
			lev.run()
			del lev

			#Move the complete telbams out of the temp_dir to the working dir
			#Only do this if we custom generated the file locations.
			if keep_in_temp:
				for source,subset_paths in output_paths.items():
					for subset,path in subset_paths.items():
						final_files.append(path)
			else:
				final_files.extend(self.__move_output_files__(output_paths))
			
			gc.collect()

		return final_files

	def __report_file_names__(self,output_paths):
		print "[Status] This run will output the following files:"
		for src,subset_paths in output_paths.items():
			for subset,output_path in subset_paths.items():
				print "\t%s" % (output_path.split("/")[-1],)
		print ""

	def __move_output_files__(self,output_paths):
		final_files = []
		for src, subset_paths in output_paths.items():
			for subset,output_path in subset_paths.items():
				try:
					move_location = output_path.replace(self._temp_dir,".")
					shutil.move(output_path,move_location) #./ being the current working dir
					final_files.append(move_location) 
				except shutil.Error,e:
					alt_filnm = "./%s_%s_%d.bam" % (src,subset,time.time()) 
					print "[Warning] Output file may already exist, you may not" \
					"have correct permissions for this file"
					print "[Update]Trying to create output using unique filename:"
					print "\t\t%s" % (alt_filnm,)
					shutil.move(output_path,alt_filnm)
					final_files.append(alt_filnm)
		return final_files

	def get_parser(self):
		#argparse imported in ./interface/parabam
		parser = self.default_parser()

		parser.add_argument('-r','--region',type=str,metavar="REGION",nargs='?',default=None
			,help="The subset process will be run only on reads from this region\n"\
			"Regions should be colon seperated as specified by samtools (eg \'chr1:1000,5000\')")
		parser.add_argument('--output','-o',metavar='OUTPUT', nargs='+',required=False
			,help="The name of the output that we wish to create. Must be same amount of space"\
			" separated entries as INPUT.")
		parser.add_argument('-s',type=int,metavar="INT",nargs='?',default=2
			,help="Further parralise subset by running this many samples side-by-side. [Default 2]")
		parser.add_argument('--debug',action="store_true",default=False,
			help="Only the first 5million reads will be processed")
		parser.add_argument('-m',action="store_true",default=False
			,help="A pair processor is used instead of a conventional processor")
		parser.add_argument('-v', choices=[0,1,2],default=0,type=int,
			help="Indicate the amount of information output by the program:\n"\
			"\t0: No output [Default]\n"\
			"\t1: Total Reads Processsed\n"\
			"\t2: Detailed output")

		return parser 

#...happily ever after