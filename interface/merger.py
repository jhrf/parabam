#!/usr/bin/env python
#Once upon a time...

#Create a bam file with all the putatative telomere reads

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
	def __init__(self,name,results,merge_type,source,destroy,total):
		super(ResultsMerge,self).__init__(name,results,destroy,total)
		self.merge_type = merge_type
		self.source = source

class HandlerMerge(parabam.Handler):

	def __init__(self,inqu,const):
		super(HandlerMerge,self).__init__(inqu,const,report=False)
		self._total = Counter()
		self._sources = const.sources
		self._merge_types = const.merge_types
		self._out_file_objects = self.__get_out_file_objects__()

	def __get_out_file_objects__(self):
		file_objects = {}
		out_file_paths = self.const.outFiles
		for src in self._sources:
			master = pysam.Samfile(self.const.master_file_path[src],"rb")
			file_objects[src] = {}
			for mrg in self._merge_types:
				cur_file_obj = pysam.Samfile(out_file_paths[src][mrg],"wb",template=master)
				file_objects[src][mrg] = cur_file_obj
			master.close()
		return file_objects

	def periodicAction(self,iterations):
		pass

	def newResultAction(self,newResult,**kwargs):
		#Handle the result. Result will always be of type framework.Result
		merge_type = newResult.merge_type
		source = newResult.source
		self._total[merge_type] = newResult.curproc #hack to record size of parent BAM

		if len(newResult.results) > 0: #Check that there are results to merge
			for result_path in newResult.results:
				result_obj = pysam.Samfile(result_path,"rb")
				for alig in result_obj.fetch(until_eof=True):
					self._out_file_objects[source][merge_type].write(alig)
				result_obj.close()
				os.remove(result_path)				

	def handlerExitFunc(self,**kwargs):
		#When the handler finishes, add the total number of reads proc'd
		#to the telbamer header. Needed for later analysis

		# DON'T KNOW IF I STILL NEED TO DO THIS.....

		# for src in self.const.sources:
		# 	curFilnm = self.const.outFiles[src][self.mrgTyp]
		# 	hedFilnm = ut.get_unique_temp("hed",tempDir=self.const.tempDir)
		# 	catFilnm = ut.get_unique_temp("cat",tempDir=self.const.tempDir)

		# 	curTelbam = pysam.Samfile(curFilnm,"rb")
		# 	newheader = curTelbam.header

		# 	totalInfo = "parent_bam_total_count:%d" % (self._total[src],)

		# 	if 'CO' in newheader:
		# 		newheader['CO'].append(totalInfo)
		# 	else:
		# 		newheader['CO'] = [totalInfo]

		# 	newTelbam = pysam.Samfile(hedFilnm,"wb",header=newheader)
		# 	curTelbam.close()
		# 	newTelbam.close()

		# 	pysam.cat("-o",catFilnm,hedFilnm,curFilnm)

		# 	os.remove(curFilnm)
		# 	os.remove(hedFilnm)
		# 	os.rename(catFilnm,curFilnm)

		for src,file_dict in self._out_file_objects.items():
			for typ,file_obj in file_dict.items():
				file_obj.close()

	#Diagnostics function
	def __storeDiag__(self,newResult):
		dirDiag = os.listdir(self.const.tempDir)

		dirDiag = [x for x in dirDiag if newResult.mrgTyp in x]
		storDiag = list(newResult.results)

		dirDiag.insert(0,"dir")
		storDiag.insert(0,"store")

		lenDiff = len(storDiag) - len(dirDiag)

		print lenDiff

		if lenDiff < 0:
			storDiag.extend(["LOL"] * abs(lenDiff))
		elif lenDiff > 0:
			dirDiag.extend(["LOL"] * lenDiff)
		else:
			"zero"

		import pprint as pp
		pp.pprint(izip(dirDiag,storDiag))