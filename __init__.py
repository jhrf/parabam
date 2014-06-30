#__init__.py

import time
import sys
import Queue as Queue2

import blessings
import pysam

from multiprocessing import Queue,Process
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

		self._updateFunc = self.__blessingsOutput__ 

	def listen(self,update):
		destroy = False
		destroyCount = 0
		dat = {} 

		iterations = 0
		stTime = time.time()
		dealt  = 0
		curproc = 0

		updateFunc = self._updateFunc #speedup alias

		if self._verbose: print "\n\n"

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

	def __blessingsOutput__(self,outstr):
		term = blessings.Terminal()
		with term.location(0,term.height-3):
			print outstr
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
		return int(time.time() - since)

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

class Processor(object):

	__metaclass__ = ABCMeta

	def __init__(self,outqu,const,debug = False):
		
		self._debug = debug
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
		collection = []
		colCount = 0

		#function alias, for optimisation
		addToCollection = self.addToCollection
		sendProc = self.__sendProc__

		bam_iterator = self.__getNextAlig__(masterBam) if not self._debug \
							else self.__getNextAligDebug__(masterBam)

		for i,alig in enumerate(bam_iterator):

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

	def __getNextAligDebug__(self,masterBam):
		for i,alig in enumerate(masterBam.fetch(until_eof=True)):
			if i < 1000000:
				yield alig
			else:
				return

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

class Leviathon(object):
	#Leviathon takes objects of processors and handlers and
	#chains them together.
	def __init__(self,processors,handlers,update=10000000):
		self._processors = processors 
		self._handlers	= handlers
		self._update = update
	
	def run(self):
		procs = []

		for p in self._processors:
			ppr = Process(target=p.run,args=(self._update,))
			ppr.start()
			procs.append(ppr)

		for h in self._handlers:
			hpr = Process(target=h.listen,args=(self._update,))
			hpr.start()
			procs.append(hpr)

		procs[-1].join() #Wait on the last handler that we started.

		map(lambda pr :pr.terminate(),procs)
		del self._processors
		del self._handlers
		del procs

class Interface(object):

	__metaclass__ = ABCMeta
	def __init__(self,tempDir):
		self._tempDir = tempDir

	@abstractmethod
	def run_cmd(self,parser,tempDir):
		#This is usualy just a function that
		#takes an argparse parser and turns 
		#passes the functions to the run function
		pass

	@abstractmethod
	def run(self):
		pass

class Const(object):
	
	def __init__(self,outFiles,tempDir,masterFnm,thresh,verbose,chunk,proc,**kwargs):
		self.outFiles = outFiles
		self.tempDir = tempDir
		self.masterFnm = masterFnm
		self.thresh = thresh
		self.verbose = verbose
		self.chunk = chunk
		self.proc = proc

		for key, val in kwargs.items():
			setattr(self,key,val)

	def addConst(**kwargs):
		for key, val in kwargs.items():
			setattr(self,key,val)

class Results(object):
	def __init__(self,name,results,destroy,curproc):
		self.name = name
		self.results = results
		self.destroy = destroy
		self.curproc = curproc

#And they all lived happily ever after...