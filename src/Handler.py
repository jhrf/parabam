#!/usr/bin/env python
#Once upon a time...

import time
import sys

import Queue as Queue2

from multiprocessing import Queue
from abc import ABCMeta, abstractmethod

class Handler(object):

	__metaclass__ = ABCMeta

	def __init__(self,inqu,const,report=True,maxDestroy=1):
		self._verbose = const.verbose
		self._filenameOut = const.outFiles
		self._inqu = inqu
		self._maxDestroy = maxDestroy
		self._report = report
		self.const = const

		self._stats = {}

		self._updateFunc = self.__decideOutput__()

	def listen(self,update):
		destroy = False
		destroyCount = 0
		dat = {} 

		iterations = 0
		stTime = time.time()
		dealt  = 0
		curproc = 0

		updateFunc = self._updateFunc #speedup alias

		while not destroy:
			iterations += 1
			#Listen for a process coming in...
			try:
				newResult = self._inqu.get(False)
				if newResult.destroy:
					destroyCount += 1
					if destroyCount == self._maxDestroy:
						destroy = True

				curproc = newResult.curproc
				if not newResult.results == {}: #Ignore empty results
					self.newResultAction(newResult)
					dealt += 1

			except Queue2.Empty:
				#Queue empty. Continue with loop
				pass

			if iterations % 100000 == 0: self.periodicAction(iterations)

			if self._verbose and self._report and iterations % update == 0:
				outstr = self.__formUpdateStr__(curproc,stTime)
				updateFunc(outstr)
				
		if self._verbose: print "[Update] All reads processed succesfully."
		self.periodicAction(iterations)
		self.handlerExitFunc()

	def __decideOutput__(self):
		if self._report and self._verbose:
			try:
				globals()["blessings"] = __import__("blessings")
				print "\n"
				return self.__blessingsOutput__
			except ImportError:
				print "[Warning] merecat will use ugly output. Please install the package ``blessings`` to correct this."
				return self.__standardOutput__
		else:
			return self.__standardOutput__

	def __blessingsOutput__(self,outstr):
		term = blessings.Terminal()
		with term.location(0,term.height-3):
			print outstr
			sys.stdout.flush()

	def __standardOutput__(self,outstr):
		if self._verbose:
			sys.stdout.write(outStr)
			sys.stdout.flush()

	def __formUpdateStr__(self,curproc,stTime):
		stats = []
		for k in self._stats:
			if type(self._stats[k]) is dict: #Account for divided stats.
				kStr = "[%s] " % (k,)
				for st in self._stats[k]:
					kStr += "%s: %d " % (st,self._stats[k][st])
				kStr = kStr[:-1] + " "

				stats.append(kStr)
			else:
				stats.append("%s: %d" % (k,self._stats[k]))

		statstr = " ||| ".join(stats)
		return "\r%s ||| [System] Active Procs: %d | Time: %d " %\
				(statstr,curproc,self.__secSince__(stTime))

	def __autoHandle__(self,res,deep=None):
		mirror = self._stats 

		if deep not in self._stats:
			self._stats[deep] = {}

		if deep:
			mirror = self._stats[deep]

		for k in res:
			if type(res[k]) is int:
				if k in mirror:
					mirror[k] += res[k]
				else:
					mirror[k] = res[k]
			elif type(res[k]) is dict:
				for k1 in res[k]:
					if type(res[k][k1]) is int:
						if k1 in mirror:
							mirror[k1] += res[k][k1]
						else:
							mirror[k1] = res[k][k1]

	def __secSince__(self,since):
		return int(time.time() - srt)

	@abstractmethod
	def periodicAction(self,iterations):
		pass

	@abstractmethod
	def newResultAction(self,newResult,**kwargs):
		#Handle the result. Result will always be of type framework.Result
		#New results action should as call autoHandle
		#self.__autoHandle__(newResult.results)
		pass

	@abstractmethod
	def handlerExitFunc(self,**kwargs):
		#When the handler finishes, do this
		pass

#...and they all lived happily ever after