#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os,argparse
import pysam
import time
import sh
import parabam

import time

import Queue as Queue2
import util as ut

from collections import Counter
from multiprocessing import Queue,Process
from itertools import izip

class ResultsMerge(parabam.Results):
	def __init__(self,name,results,subset_type,source,destroy,total,time_added):
		super(ResultsMerge,self).__init__(name,results,destroy,total)
		self.subset_type = subset_type
		self.source = source
		self.time_added = time_added

class HandlerMerge(parabam.Handler):

	def __init__(self,inqu,const):
		super(HandlerMerge,self).__init__(inqu,const,report=False)
		self._total = Counter()
		self._sources = const.sources
		self._subset_types = const.subset_types
		self._out_file_objects = self.__get_out_file_objects__()
		self._merged = 0
		self._chilling = time.time()

	def __get_out_file_objects__(self):
		file_objects = {}
		out_file_paths = self.const.outFiles
		for src in self._sources:
			master = pysam.Samfile(self.const.master_file_path[src],"rb")
			file_objects[src] = {}
			for mrg in self._subset_types:
				cur_file_obj = pysam.Samfile(out_file_paths[src][mrg],"wb",template=master)
				file_objects[src][mrg] = cur_file_obj
			master.close()
		return file_objects

	def periodicAction(self,iterations):
		pass

	def newResultAction(self,newResult,**kwargs):
		#Handle the result. Result will always be of type framework.Result
		subset_type = newResult.subset_type
		source = newResult.source
		self._total[subset_type] = newResult.curproc #hack to record size of parent BAM

		if len(newResult.results) > 0: #Check that there are results to merge
			#print "--"
			#print "since added to queue: %d | since last merge operation: %d | now merging: %d" % \
			(int(time.time() - newResult.time_added),int(time.time() - self._chilling),self._merged)
			for result_path in newResult.results:
				result_obj = pysam.Samfile(result_path,"rb")
				for alig in result_obj.fetch(until_eof=True):
					self._out_file_objects[source][subset_type].write(alig)
				result_obj.close()
				os.remove(result_path)	
			self._merged += 1		
			self._chilling = time.time()

	def add_source_to_header(self):
		#This could be moved before processing to cut down on the pysam cat overhead.
		for source in self._sources:
			for subset in self._subset_types:
				master_path = self.const.master_file_path[source]
				telbam_path = self.const.outFiles[source][subset]
				
				telbam_object = pysam.Samfile(telbam_path,"rb")	
				new_header = telbam_object.header
				telbam_object.close()
				
				source_info = "parabam_source:%s" % (master_path,)
				if 'CO' in new_header:
					new_header['CO'].append(source_info)
				else:
					new_header['CO'] = [source_info]

				header_temp_path = ut.get_unique_temp("header",tempDir=self.const.tempDir)
				header_temp_object = pysam.Samfile(header_temp_path,"wb",header=new_header)
				header_temp_object.close()

				cat_temp_path = ut.get_unique_temp("cat",tempDir=self.const.tempDir)
				pysam.cat("-o",cat_temp_path,header_temp_path,telbam_path)

				os.remove(header_temp_path)
				os.remove(telbam_path)
				os.rename(cat_temp_path,telbam_path)

	def close_all_out_files(self):
		for source,file_dict in self._out_file_objects.items():
			for subset,file_obj in file_dict.items():
				file_obj.close()
		del self._out_file_objects

	def handlerExitFunc(self,**kwargs):
		self.close_all_out_files()
		self.add_source_to_header()

