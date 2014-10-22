#__init__.py

import time
import sys
import Queue as Queue2
import gc

import pysam

import multiprocessing
from abc import ABCMeta, abstractmethod
import resource

class Task(multiprocessing.Process):

	__metaclass__ = ABCMeta

	def __init__(self,taskSet,outqu,curproc,destroy,const):
		multiprocessing.Process.__init__(self)
		self._taskSet = taskSet
		self._outqu = outqu
		self._curproc = curproc
		self._tempDir = const.tempDir
		self._destroy = destroy
		self.const = const

	def run(self):
		start_useage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		results = {}
		if len(self._taskSet) > 0: #Empty tasks sent at destory sometimes
			results = self.produceResultsDict()
			results["total"] = len(self._taskSet)

		self._outqu.put(Results(name=self.name,
								results=results,
								destroy=self._destroy,
								curproc=self._curproc))

		fin_useage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

		#Uncomment below for memory useage stats
		#print "Proc Useage %s (start/finish) : %d / %d" % (self.pid,start_useage,fin_useage)

	def __mkTmpName__(self,typ):
		#self.pid ensures that the temp names are unique.
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

		#self._term = blessings.Terminal()
		self._updateFunc = self.__standarOutput__ 

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
				
		if self._verbose and self._report: updateFunc("\n[Update] All reads processed succesfully.\n")
		self.periodicAction(iterations)
		self.handlerExitFunc()

	def __standarOutput__(self,outstr):
		if self._verbose:
			sys.stdout.write("\r" + outstr)
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
		return "\r%s || [System] Active Procs: %d | Time: %d " %\
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

	def periodicAction(self,iterations):
		#Overwrite with something that needs
		#to be done occasionally. For example
		#merging telbams.
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

	def __init__(self,outqu,const,TaskClass,task_args,debug = False):
		
		self._debug = debug
		self._master_file_path = const.master_file_path
		self._verbose = const.verbose
		self._outqu = outqu
		self._chunk = const.chunk
		self._thresh = const.thresh
		self._proc = const.proc
		self._tempDir = const.tempDir
		self.const = const
		#class which we wish to run on each read
		self._TaskClass = TaskClass
		#any additional arguments we wish to pass to the task
		self._task_args = task_args 
		self._activeProcs = []

	#Find data pertaining to assocd and all reads 
	#and divide pertaining to the chromosome that it is aligned to
	def run(self,update):
	
		self.preProcActivity(self._master_file_path)
		masterBam = self.__getMasterBam__(self._master_file_path)
		collection = []
		colCount = 0

		#function alias, for optimisation
		addToCollection = self.addToCollection
		sendProc = self.__sendProc__

		#Insert a debug iterator into the processor, this workaround doesn't
		#slow processing when not in debug mode.
		bam_iterator = self.__getNextAlig__(masterBam) if not self._debug \
							else self.__getNextAligDebug__(masterBam)

		collect_time = time.time()

		for i,alig in enumerate(bam_iterator):

			addToCollection(masterBam,alig,collection)
			colCount += 1

			if colCount == self._chunk:
				'''print "+"
				print "+ Time: %d Name#1: %s Name#2 %s Name#N %s" % (int(time.time() - collect_time),
																	collection[0].qname,
																	collection[1].qname,
																	collection[-1].qname)'''
				self.__waitOnActiveProc__(self._activeProcs,self._proc-2) # -2 for proc and handler 
				sendProc(collection)					 				   # already running
				collection = []
				colCount = 0

				collect_time = time.time()

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
		terminated_procs = []
		for pr in activeProcs:
			if not pr.is_alive():
				pr.terminate()
				terminated_procs.append(pr)
		
		for pr in terminated_procs:
			activeProcs.remove(pr)

		del terminated_procs
		#Force a collection, hopefully free memory from processes
		gc.collect() 
				
	def __sendProc__(self,collection,destroy=False):
		args = [collection,self._outqu,len(self._activeProcs)+1,destroy,self.const]
		args.extend(self._task_args)
		task = self._TaskClass(*args)
		task.start()
		self._activeProcs.append(task)

	def __getMasterBam__(self,master_file_path):
		return pysam.Samfile(master_file_path,"rb") #Open telbam for analysis

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

		map(lambda pr :pr.terminate(),procs)
		del self._processors
		del self._handlers
		del procs

class Interface(object):

	__metaclass__ = ABCMeta
	def __init__(self,tempDir,exe_dir):
		self._tempDir = tempDir
		self._exe_dir = exe_dir

	def __reportFileNames__(self,outFiles):
		print "[Update] This run will output the following files:"
		for src,mrgTypDict in outFiles.items():
			for mrgTyp,outFilePath in mrgTypDict.items():
				print "\t%s" % (outFilePath.split("/")[-1],)
		print ""

	def __moveOutputFiles__(self,outFiles):
		for src, mrgTypDict in outFiles.items():
			for mrgTyp,outFilePath in mrgTypDict.items():
				try:
					shutil.move(outFilePath,"./") #./ being the current working dir
				except shutil.Error,e:
					alt_filnm = "./%s_%s_%d.bam" % (src,mrgTyp,time.time()) 
					print "[Error] Output file may already exist, you may not" \
					"have correct permissions for this file"
					print "[Status]Trying to create using unique filename:"
					print "\t\t%s" % (alt_filnm,)
					shutil.move(outFilePath,alt_filnm)

	def __getGroup__(self,bams,names,multi):
		for i in xrange(0,len(bams),multi):
			yield (bams[i:i+multi],names[i:i+multi])

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
	
	def __init__(self,outFiles,tempDir,master_file_path,verbose,chunk,proc,**kwargs):
		self.outFiles = outFiles
		self.tempDir = tempDir
		self.master_file_path = master_file_path
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