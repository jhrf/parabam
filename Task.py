#!/usr/bin/env python
#Once upon a time...

from multiprocessing import Process
from abc import ABCMeta, abstractmethod

class Task(Process):

	__metaclass__ = ABCMeta

	def __init__(self,taskSet,outqu,curproc,destroy,const):
		Process.__init__(self)

		self._taskSet = taskSet
		self._outqu = outqu
		self._curproc = curproc
		self._tempDir = const.tempDir
		self._destroy = destroy
		self.const = const

	def run(self):
		results = {}
		if len(self._taskSet) > 0: #Empty tasks sent at destory sometimes
			results = self.produceResultsDict()
			results["total"] = len(self._taskSet)

		self._outqu.put(Results(name=self.name,
								results=results,
								destroy=self._destroy,
								curproc=self._curproc))

	def __mkTmpName__(self,typ):
		#---FIX---
		#If run in parallel this could create two files of 
		#exactly the same name, although unlikely.
		return "%s/%s_%s_TmpTelbam.bam" % (self._tempDir,typ,self.pid)

	def __getMate__(self,master,alig):
		strtPos = master.tell()
		mate = master.mate(alig)
		master.seek(strtPos)
		return mate

	@abstractmethod
	def produceResultsDict(self):
		#Must return dict
		pass

#And they all lived happily ever after...