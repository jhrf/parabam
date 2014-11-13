#__init__.py

import os
import time
import datetime
import sys
import Queue as Queue2
import gc
import shutil
import argparse

import pysam

import multiprocessing
from abc import ABCMeta, abstractmethod
import resource

#Tasks are started by parabam.Processor. Once started 
#they will carryout a predefined task on the provided
#subset of reads (task_set).
class Task(multiprocessing.Process):

	__metaclass__ = ABCMeta

	def __init__(self,task_set,outqu,curproc,destroy,const):
		multiprocessing.Process.__init__(self)
		self._task_set = task_set
		self._outqu = outqu
		self._curproc = curproc
		self._temp_dir = const.temp_dir
		self._destroy = destroy
		self.const = const

	def run(self):
		results = {}
		results = self.__generate_results__()
		results["total"] = len(self._task_set)

		self._outqu.put(Package(name=self.name,
								results=results,
								destroy=self._destroy,
								curproc=self._curproc))

		#Trying to control mem useage
		del self._task_set 
		gc.collect()

	def __get_temp_path__(self,typ):
		#self.pid ensures that the temp names are unique.
		return "%s/%s_%s_parabam_temp.bam" % (self._temp_dir,typ,self.pid)

	def __get_mate__(self,master,alig):
		start = master.tell()
		mate = master.mate(alig)
		master.seek(start)
		return mate

	#Generate a dictionary of results using self._task_set
	@abstractmethod
	def __generate_results__(self):
		#Should always return a dictionary
		pass

#Handlers recieves the results of parabam.Task via
#a multiprocessing Queue.
class Handler(object):

	__metaclass__ = ABCMeta

	def __init__(self,inqu,const,report=True,destroy_limit=1):
		self._inqu = inqu
		self._destroy_limit = destroy_limit
		self._report = report
		self.const = const

		self._verbose = const.verbose
		self._output_paths = const.output_paths

		self._stats = {}

		self._update_output = self.__standard_output__ 

	def listen(self,update_interval):
		destroy = False
		destroy_count = 0

		iterations = 0
		start_time = time.time()
		dealt  = 0
		curproc = 0

		update_output = self._update_output #speedup alias
		while not destroy:
			iterations += 1
			#Listen for a process coming in...
			try:
				new_package = self._inqu.get(False)
				if new_package.destroy:
					#destroy limit allows for multiple processors
					destroy_count += 1
					if destroy_count == self._destroy_limit:
						destroy = True

				curproc = new_package.curproc 
				if not new_package.results == {}:#If results are present...
					self.__new_package_action__(new_package) #Handle the results
					dealt += 1

			except Queue2.Empty:
				#Queue empty. Continue with loop
				pass

			if iterations % 100000 == 0: 
				self.__periodic_action__(iterations)
				gc.collect()

			if self._verbose and self._report and iterations % update_interval == 0:
				outstr = self.__format_update__(curproc,start_time)
				update_output(outstr)

		self.__periodic_action__(iterations)
		if self._verbose and self._report: update_output("\n[output_Update] All reads processed succesfully.\n")
		self._inqu.close()
		self.__handler_exit__()

	def __standard_output__(self,outstr):
		if self._verbose:
			sys.stdout.write("\r" + outstr)
			sys.stdout.flush()

	def __format_update__(self,curproc,start_time):
		stats = []
		for stat_name in self._stats:
			if type(self._stats[stat_name]) is dict: #Account for divided stats.
				stat_update_str = "[%s] " % (stat_name,)
				for data in self._stats[stat_name]:
					stat_update_str += "%s: %d " % (data,self._stats[stat_name][data])
				stat_update_str = stat_update_str[:-1] + " "

				stats.append(stat_update_str)
			else:
				stats.append("%s: %d" % (stat_name,self._stats[stat_name]))

		statstr = " | ".join(stats)
		return "\r%s - [System] Active Procs: %d | Time: %d " %\
				(statstr,curproc,self.__secs_since__(start_time))

	#This code is a little ugly. Essentially, given a results
	#dictionary, it will go through and create a sensible output
	#string.
	def __auto_handle__(self,results,deep=None):
		stats = self._stats 

		if deep not in self._stats:
			self._stats[deep] = {}

		if deep:
			stats = self._stats[deep]

		for key_one in results:
			if type(results[key_one]) is int:
				if key_one in stats:
					stats[key_one] += results[key_one]
				else:
					stats[key_one] = results[key_one]
			elif type(results[key_one]) is dict:
				for key_two in results[key_one]:
					if type(results[key_one][key_two]) is int:
						if key_two in stats:
							stats[key_two] += results[key_one][key_two]
						else:
							stats[key_two] = results[key_one][key_two]

	def __secs_since__(self,since):
		return int(time.time() - since)

	def __periodic_action__(self,iterations):
		#Overwrite with something that needs
		#to be done occasionally.		
		pass

	@abstractmethod
	def __new_package_action__(self,new_package,**kwargs):
		#Handle the output of a task. Input will always be of type
		#parabam.Package. New results action should always call
		#self.__auto_handle__(new_package.results)
		pass

	@abstractmethod
	def __handler_exit__(self,**kwargs):
		#When the handler finishes, do this
		pass

#The Processor iterates over a BAM, subsets reads and
#then starts a parbam.Task process on the subsetted reads
class Processor(object):

	__metaclass__ = ABCMeta

	def __init__(self,outqu,const,TaskClass,task_args,debug=False):

		self._debug = debug
		self._outqu = outqu

		#Constants
		self._master_file_path = const.master_file_path
		self._verbose = const.verbose
		self._chunk = const.chunk
		self._thresh = const.thresh
		self._proc = const.proc
		self._temp_dir = const.temp_dir
		self.const = const

		#Class which we wish to run on each read
		self._TaskClass = TaskClass
		self._task_args = task_args 
		
		self._active_tasks = []

	#Find data pertaining to assocd and all reads 
	#and divide pertaining to the chromosome that it is aligned to
	def run(self,update):
	
		self.__pre_processor__(self._master_file_path)
		master_bam = self.__get_master_bam__(self._master_file_path)
		collection = []
		collection_count = 0

		#function alias, for optimisation
		add_to_collection = self.__add_to_collection__
		start_task = self.__start_task__

		#Insert a debug iterator into the processor, this workaround doesn't
		#slow processing when not in debug mode.
		bam_iterator = self.__get_next_alig__(master_bam) if not self._debug \
							else self.__get_next_alig_debug__(master_bam)

		for i,alig in enumerate(bam_iterator):

			add_to_collection(master_bam,alig,collection)
			collection_count += 1

			if collection_count == self._chunk:
				self.__wait_for_tasks__(self._active_tasks,self._proc) # -2 for proc and handler 
				start_task(collection)					 				   # already running
				collection = []
				collection_count = 0

		self.__wait_for_tasks__(self._active_tasks,0)
		start_task(collection,destroy=True)

		self.__end_processing__(master_bam)

	def __output__(self,outstr):
		if self._verbose:
			sys.stdout.write(outstr)
			sys.stdout.flush()

	#__wait_for_tasks__ controls the amount of currently running processors
	#it waits until a processor is free and then allows the creation of a new
	#process
	def __wait_for_tasks__(self,active_tasks,max_tasks):
		update_tasks = self.__update_tasks__ #optimising alias
		update_tasks(active_tasks)
		currently_active = len(active_tasks)

		if max_tasks > currently_active:
			return

		i = 0
		report = False

		while(max_tasks < currently_active):
			currently_active = len(active_tasks)
			update_tasks(active_tasks)
			i += 1
			if i % 100000 == 0 and self._verbose:
				if not report: #Only say the following once
					if max_tasks > 0:
						sys.stdout.write("\r[Status] Processors all busy...waiting")
						pass
					elif max_tasks == 0:
						sys.stdout.write("\r[Status] Waiting for current tasks to finish before sending last batch")
				time.sleep(5)
				report = True

		if self._verbose and report: 
			if max_tasks == 0: 
				print "\r[Update] All tasks have finished, sending final task"
			elif report:
				sys.stdout.write("\r                                                           ")
				sys.stdout.write("\r[Status] Normal functioning resumed")

	def __update_tasks__(self,active_tasks):
		terminated_procs = []
		for pr in active_tasks:
			if not pr.is_alive():
				pr.terminate()
				terminated_procs.append(pr)
		
		for pr in terminated_procs:
			active_tasks.remove(pr)

		del terminated_procs
		#Force a collection, hopefully free memory from processes
		gc.collect() 
				
	def __start_task__(self,collection,destroy=False):
		args = [collection,self._outqu,len(self._active_tasks)+1,destroy,self.const]
		args.extend(self._task_args)
		task = self._TaskClass(*args)
		task.start()
		self._active_tasks.append(task)

	def __get_master_bam__(self,master_file_path):
		return pysam.Samfile(master_file_path,"rb") #Open telbam for analysis

	#Should be over written if we don't want data from BAM file.
	def __get_next_alig__(self,master_bam):
		for alig in master_bam.fetch(until_eof=True):
			yield alig

	def __get_next_alig_debug__(self,master_bam):
		for i,alig in enumerate(master_bam.fetch(until_eof=True)):
			if i < 10000000:
				yield alig
			else:
				return

	def __end_processing__(self,master):
		master.close()
		for proc in self._active_tasks:
			proc.join()
			proc.terminate()
		del self._active_tasks
		gc.collect()

	@abstractmethod
	def __pre_processor__(self,master_bam):
		pass

	@abstractmethod
	def __add_to_collection__(self,master,item,collection):
		pass

class Leviathon(object):
	#Leviathon takes objects of processors and handlers and
	#chains them together.
	def __init__(self,processors,handlers,update=10000000):
		self._processors = processors 
		self._handlers	= handlers
		self._update = update
	
	def run(self):
		#We spawn processes rather than fork to save memory
		procs = []

		for p in self._processors:
			ppr = multiprocessing.Process(target=p.run,args=(self._update,))
			ppr.start()
			procs.append(ppr)

		for h in self._handlers:
			hpr = multiprocessing.Process(target=h.listen,args=(self._update,))
			hpr.start()
			procs.append(hpr)

		procs[-1].join() #Wait on the last handler that we started.

		for proc in procs:
			proc.join()
			proc.terminate()

		del self._processors
		del self._handlers
		del procs

#Provides a conveinant way for providing an Interface to parabam
#programs. Includes default command_args and framework for 
#command-line and programatic invocation. 
class Interface(object):

	__metaclass__ = ABCMeta
	def __init__(self,temp_dir,exe_dir):
		self._temp_dir = temp_dir
		self._exe_dir = exe_dir

	def __introduce__(self,name):
		intro =  "%s has started. Start Time: " % (name,)\
			+ datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
		underline = "-" * len(intro)
		print intro
		print underline

	def __goodbye__(self,name):
		print "%s has finished. End Time: " % (name,)\
			+ datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

	def __get_group__(self,bams,names,multi):
		for i in xrange(0,len(bams),multi):
			yield (bams[i:i+multi],names[i:i+multi])

	def __get_basename__(self,path):
		base = os.path.basename(path)
		if "." in base:
			base = base.rpartition(".")[0]
		return base

	def __sort_and_index__(self,fnm,verbose=False,tempDir=".",name=False):
		if not os.path.exists(fnm+".bai"):
			if verbose: print "%s is not indexed. Sorting and creating index file." % (fnm,)
			tempsort_path = get_unqiue_tmp_path("SORT",tempDir=tempDir)
			optStr = "-nf" if name else "-f"
			pysam.sort(optStr,fnm,tempsort_path)
			os.remove(fnm)
			os.rename(tempsort_path,fnm)
			if not name: pysam.index(fnm)

	def default_parser(self):

		parser = argparse.ArgumentParser(conflict_handler='resolve')

		parser.add_argument('-p',type=int,nargs='?',default=8
			,help="The amount of processors you wish the program to use."\
					" Ignored if argument 's' is provided")
		parser.add_argument('-c',type=int,nargs='?',default=25000
			,help="How many arguments the preprocesor should store before"\
			"sending to the analysis module")
		parser.add_argument('-v',action="store_true",default=False
			,help='If set the program will output more information to terminal')

		return parser

	@abstractmethod
	def run_cmd(self,parser):
		#This is usualy just a function that
		#takes an argparse parser and turns 
		#passes the functions to the run function
		pass

	@abstractmethod
	def run(self):
		pass

class Const(object):
	
	def __init__(self,output_paths,temp_dir,master_file_path,verbose,chunk,proc,**kwargs):
		self.output_paths = output_paths
		self.temp_dir = temp_dir
		self.master_file_path = master_file_path
		self.verbose = verbose
		self.chunk = chunk
		self.proc = proc

		for key, val in kwargs.items():
			setattr(self,key,val)

class Package(object):
	def __init__(self,name,results,destroy,curproc):
		self.name = name
		self.results = results
		self.destroy = destroy
		self.curproc = curproc

#And they all lived happily ever after...