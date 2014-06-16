#!/usr/bin/env python
#Once upon a time...

from multiprocessing import Process

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

#...and they all lived happily ever after