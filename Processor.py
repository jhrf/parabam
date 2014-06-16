#!/usr/bin/env python
#Once upon a time...

import pysam
import time
import sys

from abc import ABCMeta, abstractmethod

class Processor(object):

	__metaclass__ = ABCMeta

	def __init__(self,outqu,const):
		
		self._masterFnm = const.masterFnm
		self._verbose = const.verbose
		self._outqu = outqu
		self._chunk = const.chunk
		self._thresh = const.thresh
		self._proc = const.proc
		self._tempDir = const.tempDir
		self.const = const
		self._activeProcs = []

	#Find data pertaining to assocd and all reads 
	#and divide pertaining to the chromosome that it is aligned to
	def run(self,update):
		self.preProcActivity(self._masterFnm)

		masterBam = self.__getMasterBam__(self._masterFnm)
		getNextAlig = self.__getNextAlig__
		collection = []
		colCount = 0

		#function alias, for optimisation
		addToCollection = self.addToCollection
		sendProc = self.__sendProc__

		for i,alig in enumerate(getNextAlig(masterBam)):

			addToCollection(masterBam,alig,collection)
			colCount += 1

			if colCount == self._chunk:
				self.__waitOnActiveProc__(self._activeProcs,self._proc-2) # -2 for proc and handler 
				sendProc(collection)					 				   # already running
				collection = []
				colCount = 0

		self.__waitOnActiveProc__(self._activeProcs,0)
		sendProc(collection,destroy=True)

		self.__endProcessing__(masterBam)


	#__waitOnActiveProc__ controls the amount of currently running processors
	#it waits until a processor is free and then allows the creation of a new
	#process
	def __waitOnActiveProc__(self,activeProcs,maxProc):
		procUpdate = self.__procUpdate__ #optimising alias
		procUpdate(activeProcs)
		activeLen = len(activeProcs)

		if maxProc > activeLen:
			return

		i = 0
		stWait = time.time()
		report = False

		while( maxProc < activeLen):
			activeLen = len(activeProcs)
			procUpdate(activeProcs)
			i += 1
			if i % 100000 == 0 and self._verbose:
				if not report: #Only say the following once
					if maxProc > 0:
						sys.stdout.write("\r[Status] Processors all busy...waiting")
						pass
					elif maxProc == 0:
						sys.stdout.write("\r[Status] Waiting for current tasks to finish before sending last batch")
				time.sleep(7)
				report = True

		if self._verbose and report: 
			if maxProc == 0: 
				print "\r[Update] All tasks have finished, sending final task"
			elif report:
				sys.stdout.write("\r                                                           ")
				sys.stdout.write("\r[Status] Normal functioning resumed")

	def __procUpdate__(self,activeProcs):
		for pr in activeProcs:
			if not pr.is_alive():
				pr.terminate()
				activeProcs.remove(pr)
				del pr			

	def __sendProc__(self,collection,destroy=False):
		#self,taskSet,outqu,curproc,destroy,const
		args = [collection,self._outqu,len(self._activeProcs)+1,destroy,self.const]
		(taskClass,customargs) = self.getCustomTask()
		args.extend(customargs)
		task = taskClass(*args)
		task.start()
		self._activeProcs.append(task)

	def __getMasterBam__(self,masterFnm):
		return pysam.Samfile(masterFnm,"rb") #Open telbam for analysis

	#Should be over written if we don't want data from BAM file.
	def __getNextAlig__(self,masterBam):
		for alig in masterBam.fetch(until_eof=True):
			yield alig

	def __endProcessing__(self,master):
		master.close()

	@abstractmethod
	def preProcActivity(self,masterBam):
		pass

	@abstractmethod		
	def getCustomTask(self):
		#Must return class and dict of customargs
		pass

	@abstractmethod
	def addToCollection(self,master,item,collection):
		pass

#...and they all lived happily ever after