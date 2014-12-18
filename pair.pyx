import os
import time
import datetime
import sys
import Queue as Queue2
import gc
import shutil

import argparse
import parabam
import pysam

import multiprocessing
from abc import ABCMeta, abstractmethod
import resource

class ProcessorSubset(parabam.tools.Processor):

	def __init__(self,object outqu,object const,object TaskClass,object task_args,object debug=False):
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

	def __init__(self,object outqu,object const,object TaskClass,object task_args,object debug=False):
		super(PairProcessor,self).__init__(outqu,const,TaskClass,task_args,debug=False)
		self._loners = {}
		
		self._loner_count = 0
		master = pysam.Samfile(self._master_file_path[task_args[0]],"rb")
		self._loners_object = pysam.Samfile("%s/loners" % (self._temp_dir,),"wb",template=master)
		master.close()

	def __add_to_collection__(self,master,item,collection):
		loners = self._loners
		cdef temp_count = self._loner_count
		if not item.is_secondary:
			try:
				mate = loners[item.qname] 
				del loners[item.qname]
				temp_count -= 1
				collection.append( (item,mate,) )
			except KeyError:
				#Could implement a system where by long standing
				#unpaired reads are stored to be run at the end 
				#of the program, otherwise we risk clogging memory
				loners[item.qname] = item
				temp_count += 1

		if temp_count > 200000:
			self._loner_count = 0
			self.__stash_loners__(loners)
			del self._loners
			gc.collect()
			self._loners = {}
		self._loner_count = temp_count

	def __stash_loners__(self,loners):
		for name,read in loners.items():
			self._loners_object.write(read)

	def __end_processing__(self,master):
		super(PairProcessor,self).__end_processing__(master)
		self.__stash_loners__(self._loners)
		self._loners_object.close()