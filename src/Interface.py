#!/usr/bin/env python

#Once upon a time...

import pysam
import pdb
import time
import sys
import sh

import Queue as Queue2
import numpy as np
import util as ut

from multiprocessing import Queue,Process
from abc import ABCMeta, abstractmethod

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
