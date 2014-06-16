#!/usr/bin/env python
#Once upon a time...

class Results(object):
	def __init__(self,name,results,destroy,curproc):
		self.name = name
		self.results = results
		self.destroy = destroy
		self.curproc = curproc