#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os,argparse
import pysam
import time
import parabam

import time

import Queue as Queue2

from collections import Counter
from multiprocessing import Queue,Process
from itertools import izip

class MergePackage(parabam.Package):
	def __init__(self,name,results,subset_type,source,destroy,total,time_added):
		super(MergePackage,self).__init__(name,results,destroy,total)
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
		output_paths = self.const.output_paths
		for src in self._sources:
			master = pysam.Samfile(self.const.master_file_path[src],"rb")
			file_objects[src] = {}
			for mrg in self._subset_types:
				cur_file_obj = pysam.Samfile(output_paths[src][mrg],"wb",template=master)
				file_objects[src][mrg] = cur_file_obj
			master.close()
		return file_objects

	def __periodic_action__(self,iterations):
		pass

	def __new_package_action__(self,new_package,**kwargs):
		#Handle the result. Result will always be of type framework.Result
		subset_type = new_package.subset_type
		source = new_package.source
		self._total[subset_type] = new_package.curproc #hack to record size of parent BAM

		if len(new_package.results) > 0: #Check that there are results to merge
			(int(time.time() - new_package.time_added),int(time.time() - self._chilling),self._merged)
			for result_path in new_package.results:
				result_obj = pysam.Samfile(result_path,"rb")
				for alig in result_obj.fetch(until_eof=True):
					self._out_file_objects[source][subset_type].write(alig)
				result_obj.close()
				os.remove(result_path)	
			self._merged += 1		
			self._chilling = time.time()

	def __add_source_to_header__(self):
		#This could be moved before processing to cut down on the pysam cat overhead.
		for source in self._sources:
			for subset in self._subset_types:
				master_path = self.const.master_file_path[source]
				telbam_path = self.const.output_paths[source][subset]
				
				telbam_object = pysam.Samfile(telbam_path,"rb")	
				new_header = telbam_object.header
				telbam_object.close()
				
				source_info = "parabam_source:%s" % (master_path,)
				if 'CO' in new_header:
					new_header['CO'].append(source_info)
				else:
					new_header['CO'] = [source_info]

				header_temp_path = self.__get_unique_temp__("header",temp_dir=self.const.temp_dir)
				header_temp_object = pysam.Samfile(header_temp_path,"wb",header=new_header)
				header_temp_object.close()

				cat_temp_path = self.__get_unique_temp__("cat",temp_dir=self.const.temp_dir)
				pysam.cat("-o",cat_temp_path,header_temp_path,telbam_path)

				os.remove(header_temp_path)
				os.remove(telbam_path)
				os.rename(cat_temp_path,telbam_path)

	def __get_unique_temp__(self,temp_type,temp_dir="."):
		return "%s/%sTEMP%d.bam" % (temp_dir,temp_type,int(time.time()),)

	def __close_all_out_files__(self):
		for source,file_dict in self._out_file_objects.items():
			for subset,file_obj in file_dict.items():
				file_obj.close()
		del self._out_file_objects

	def __handler_exit__(self,**kwargs):
		self.__close_all_out_files__()
		self.__add_source_to_header__()

#...happily ever after