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

cdef class Handler:

	cdef int _destroy_limit

	def __init__(self,object inqu,object const,object report=True,int destroy_limit=1):
		self._inqu = inqu
		self._destroy_limit = destroy_limit
		self._report = report
		self.const = const

		self._verbose = const.verbose
		self._output_paths = const.output_paths

		self._stats = {}

		if const.verbose == 1:
			self._verbose = True
			self._update_output = self.__level_1_output__
		elif const.verbose == 2 or const.verbose == True:
			#In the case verbose is simply "True" or "level 2"
			self._verbose = True
			self._update_output = self.__standard_output__
		else:#catching False and -v0
			self._verbose = False
			self._update_output = self.__standard_output__

	def __level_1_output__(self,out_str):
		total_procd = self.__get_total_processed_reads__()
		time = out_str.partition("Time: ")[2]
		sys.stdout.write("[Update] Processed: %d Time: %s\n" % (total_procd,time))
		sys.stdout.flush()

	def __standard_output__(self,outstr):
		sys.stdout.write("\r" + outstr)
		sys.stdout.flush()

	def listen(self,update_interval):
		destroy = False
		cdef int destroy_count = 0

		cdef int iterations = 0
		cdef int start_time = time.time()
		cdef int dealt  = 0
		cdef int curproc = 0

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
				time.sleep(1)

			if iterations % 10 == 0: 
				self.__periodic_action__(iterations)

			if self._verbose and self._report and iterations % update_interval == 0:
				outstr = self.__format_update__(curproc,start_time)
				update_output(outstr)

		self.__periodic_action__(iterations)
		self._inqu.close()
		self.__handler_exit__()

	def __format_update__(self,curproc,start_time):
		stats = []

		for stat_name in self._stats:
			if type(self._stats[stat_name]) is dict: #Account for divided stats.
				stat_update_str = "[%s] " % (stat_name,)
				for data in self._stats[stat_name]:
					stat_update_str += "%s:%d " % (data,self._stats[stat_name][data])
				stat_update_str = stat_update_str[:-1] + " "

				stats.append(stat_update_str)
			else:
				stats.append("%s: %d" % (stat_name,self._stats[stat_name]))

		statstr = " ".join(stats)
		return "\r%s | Tasks: %d Time: %d " %\
				(statstr,curproc,self.__secs_since__(start_time),)

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

	def __new_package_action__(self,new_package,**kwargs):
		#Handle the output of a task. Input will always be of type
		#parabam.Package. New results action should always call
		#self.__auto_handle__(new_package.results)
		pass

	def __handler_exit__(self,**kwargs):
		#When the handler finishes, do this
		pass

#The Processor iterates over a BAM, subsets reads and
#then starts a parbam.Task process on the subsetted reads
cdef class Processor:

	cdef int _chunk,_proc,_verbose

	def __init__(self,object outqu,object const,object TaskClass,list task_args,object debug=False):

		self._debug = debug

		self._chunk = const.chunk
		self._proc = const.proc

		self._master_file_path = const.master_file_path
		self._verbose = const.verbose
		self._temp_dir = const.temp_dir

		self._outqu = outqu
		self.const = const

		#Class which we wish to run on each read
		self._TaskClass = TaskClass
		self._task_args = task_args 

		self._active_tasks = []

	#Find data pertaining to assocd and all reads 
	#and divide pertaining to the chromosome that it is aligned to
	def run(self,int update):
	
		self.__pre_processor__(self._master_file_path)
		
		cdef object master_bam = self.__get_master_bam__(self._master_file_path)
		cdef list collection = []
		cdef int collection_count = 0

		#function alias, for optimisation
		add_to_collection = self.__add_to_collection__
		start_task = self.__start_task__
		wait_for_tasks = self.__wait_for_tasks__

		#Insert a debug iterator into the processor, this workaround doesn't
		#slow processing when not in debug mode.
		bam_iterator = self.__get_next_alig__(master_bam) if not self._debug \
							else self.__get_next_alig_debug__(master_bam)

		for i,alig in enumerate(bam_iterator):

			add_to_collection(master_bam,alig,collection)
			collection_count += 1

			if collection_count == self._chunk:
				wait_for_tasks(self._active_tasks,self._proc) # -2 for proc and handler 
				start_task(collection)		
				del collection
				collection = []			 				   # already running
				collection_count = 0

		wait_for_tasks(self._active_tasks,0)
		start_task(collection,destroy=True)

		self.__end_processing__(master_bam)

	def __output__(self,outstr):
		if self._verbose:
			sys.stdout.write(outstr)
			sys.stdout.flush()

	#__wait_for_tasks__ controls the amount of currently running processors
	#it waits until a processor is free and then allows the creation of a new
	#process
	def __wait_for_tasks__(self,list active_tasks,int max_tasks):
		update_tasks = self.__update_tasks__ #optimising alias
		update_tasks(active_tasks)
		cdef int currently_active = len(active_tasks)

		if max_tasks > currently_active:
			return

		while(max_tasks < currently_active):
			update_tasks(active_tasks)
			currently_active = len(active_tasks)
			time.sleep(1)

		return 

	def __update_tasks__(self,active_tasks):
		terminated_procs = []
		for pr in active_tasks:
			if not pr.is_alive():
				pr.terminate()
				terminated_procs.append(pr)
		
		for pr in terminated_procs:
			active_tasks.remove(pr)

		if len(terminated_procs) > 0:
			#Only invoked the GC if there is the possibility of
			#collecting any memory.
			del terminated_procs
			gc.collect()
					
	def __start_task__(self,collection,destroy=False):
		args = [collection,self._outqu,len(self._active_tasks),destroy,self.const]
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
			if i < 5000000:
				yield alig
			else:
				return

	def __end_processing__(self,master):
		master.close()
		for proc in self._active_tasks:
			proc.join()
			proc.terminate()
		#del self._active_tasks
		#gc.collect()

	def __pre_processor__(self,master_bam):
		pass

	def __add_to_collection__(self,master,item,collection):
		pass
