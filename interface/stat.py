import pysam
import parabam
import shutil
import time
import sys
import os
import copy
import numpy as np

from multiprocessing import Queue

from parabam.interface import util as ut
from parabam.interface import merger

from abc import ABCMeta, abstractmethod

class TaskStat(parabam.Task):

	__metaclass__ = ABCMeta

	def __init__(self,taskSet,outqu,curproc,destroy,const,source):
		super(TaskStat, self).__init__(taskSet=taskSet,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
		self._counts = {}
		self._local_structures = {}

	def produceResultsDict(self):
		master = pysam.Samfile(self._master_file_path,"rb")

		for name,struct in self.const.user_structures.items():
			self._local_structures[name] = struct.empty_clone()
			self._counts[name] = 0

		self.handle_task_set(self._taskSet,master)

		master.close() #Close the master bamfile

		results = {}
		results["source"] = self._source
		results["structures"] = self.unpack_structures(self._local_structures)
		results["counts"] = self._counts

		return results

	def unpack_structures(self,structures):
		unpacked = []
		for name,struc in structures.items():
			unpacked.append( (name,struc.data) )

		return unpacked

	def handle_task_set(self,task_set,master):
		engine = self.const.user_engine
		constants = self.const.user_constants
		local_structures = self._local_structures
		for alig in task_set:
			results = engine(alig,constants,master)
			if results:
				for name,package in results.items():
					self._counts[name] += 1
					local_structures[name].add(**package)

class HandlerStat(parabam.Handler):

	def __init__(self,inqu,const):
		super(HandlerStat,self).__init__(inqu,const,maxDestroy=len(const.sources))

		self._sources = const.sources
		self._final_structures = {}

		for source in self._sources:
			self._final_structures[source] = {}
			for name,struc in self.const.user_structures.items():
				self._final_structures[source][struc.name] = struc.empty_clone() 

	def newResultAction(self,newResult,**kwargs):
		results = newResult.results
		source = results["source"]
		self.__autoHandle__(results,source)
		for name,data in results["structures"]:
			final_struc = self._final_structures[source][name]
			final_struc.merge(data)

	def periodicAction(self,iterations):
		pass

	def handlerExitFunc(self,**kwargs):
		for source in self._sources:
			for name,struc in self._final_structures[source].items():
				struc.write_to_csv(self.const.outFiles,source,self.const.outmode)

class ProcessorStat(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args):
		super(ProcessorStat,self).__init__(outqu,const,TaskClass,task_args)
		self._source = task_args[0] #Defined in the run function within Interface

	def __getMasterBam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def addToCollection(self,master,alig,collection):
		collection.append(alig)

	def preProcActivity(self,master_file_path):
		pass

class UserStructure(object):

	__metaclass__ = ABCMeta

	def __init__(self,name,struc_type,store_method,data):
		self.struc_type = struc_type
		self.store_method = store_method 
		self.data = data
		self.org_data = copy.copy(data)
		self.name = name

		if store_method == "max":
			self.add = self.add_max
			self.merge = self.merge_max
		elif store_method == "min":
			self.add = self.add_min
			self.merge = self.merge_min
		else:
			self.add = self.add_cumu
			self.merge = self.merge_cumu

	def max_decision(self,result,existing):
		return max(result,existing)

	def min_decision(self,result,exisiting):
		return min(result,exisiting)

	@abstractmethod
	def empty_clone(self):
		pass

	@abstractmethod
	def add_max(self,result):
		pass

	@abstractmethod
	def add_min(self,result):
		pass

	@abstractmethod
	def add_cumu(self,result):
		pass

	@abstractmethod
	def merge_max(self,result):
		pass

	@abstractmethod
	def merge_min(self,result):
		pass

	@abstractmethod
	def merge_cumu(self,result):
		pass

	@abstractmethod
	def write_to_csv(self,source):
		pass

class NumericStructure(UserStructure):
	def __init__(self,name,struc_type,store_method,data,log_scaling=False):
		super(NumericStructure,self).__init__(name,struc_type,store_method,data)
		self.log_scaling = log_scaling

	def empty_clone(self):
		return NumericStructure(self.name,self.struc_type,self.store_method,self.org_data)

	def add_cumu(self,result):
		self.data += result

	def add_max(self,result):
		self.data = self.max_decision(result,self.data)

	def add_min(self,result):
		self.data = self.min_decision(result,self.data)

	def merge_max(self,result):
		self.data = self.max_decision(self.data,result)
		del result

	def merge_min(self,result):
		self.data = self.min_decision(self.data,result)
		del result

	def merge_cumu(self,result):
		self.data += result
		del result

	def write_to_csv(self,out_paths,source,mode):
		if mode == "a":
			first_col=source
		elif mode == "s":
			first_col=self.name
		write_line = "%s,%.3f\n" % (first_col,self.data,)

		with open(out_paths[source][self.name],"a") as out_object:
			out_object.write(write_line)

class ArrayStructure(UserStructure):
	def __init__(self,name,struc_type,store_method,data):
		super(NumericStructure,self).__init__(name,struc_type,store_method,data)
		self.log_scaling = log_scaling

	def empty_clone(self):
		return NumericStructure(self.name,self.struc_type,self.store_method,self.org_data)

	def add_max(self,result,coords):
		existing = data[coords]
		self.data[coords] = self.max_decision(result,existing)

	def add_min(self,result,coords):
		existing = data[coords]
		self.data[coords] = self.min_decision(result,existing)

	def add_cumu(self,result,coords):
		self.data[coords] += result

	def merge_max(self,result):
		self.data = np.maximum(self.data,result)
		del result

	def merge_min(self,result):
		self.data = np.minimum(self.data,result)
		del result

	def merge_cumu(self,result):
		self.data += result
		del result

	def write_to_csv(self,out_paths,source,mode):
		np.savetxt(out_paths[source][self.name],self.data,fmt="%.3f",delimiter=",")

class Interface(parabam.Interface):

	def __init__(self,tempDir,exe_dir):
		super(Interface,self).__init__(tempDir,exe_dir)
	
	def run_cmd(self,parser):

		cmd_args = parser.parse_args()

		verbose = cmd_args.v
		module = __import__(cmd_args.instruc, fromlist=[''])
		user_engine = module.engine
		user_constants = {}
		user_struc_blueprint = {}
		module.set_constants(user_constants)
		module.set_structures(user_structures_blueprint)

		self.run(
			input_bams=cmd_args.input,
			proc= cmd_args.p,
			chunk= cmd_args.c,
			verbose= verbose,
			user_constants = user_constants,
			user_engine = user_engine,
			user_struc_blueprint = user_struc_blueprint,
			engine_is_class = False,
			outmode = cmd_args.outmode)
	
	def run(self,input_bams,proc,chunk,verbose,user_constants,
			user_engine,engine_is_class,user_struc_blueprint,outmode,multi=4):

		user_structures = self.create_structures(user_struc_blueprint)
		outputs = [ ut.get_bam_basename(b) for b in input_bams ]
		for input_group,output_group in self.__getGroup__(input_bams,outputs,multi=multi):
			
			qu = Queue()

			master_file_path = {}

			for master,source in zip(input_group,output_group):
				master_file_path[source] = master
					
			procrs = []
			handls = []

			out_paths = self.create_out_paths(outmode,output_group,user_structures)

			const = parabam.Const(outFiles=out_paths,tempDir=self._tempDir,
								master_file_path=master_file_path,
								chunk=chunk,proc=(proc // len(input_group)),
								verbose=verbose,thresh=0,
								sources=output_group,
								exe_dir=self._exe_dir,
								user_constants=user_constants,
								user_structures=user_structures,
								user_engine=user_engine,
								outmode=outmode)

			for source in output_group:
				if engine_is_class:
					if not issubclass(user_engine,TaskStat):
						raise Exception("[ERROR]\tThe class provided to parabam must be a subclass of\n"\
										"\tparabam.interface.subset.TaskStat. Please consult the parabam manual.")
					cur_args = [source]
					procrs.append(ProcessorStat(outqu=qu,
											const=const,
											TaskClass=user_engine,
											task_args=cur_args))
				else:
					procrs.append(ProcessorStat(outqu=qu,
												const=const,
												TaskClass=TaskStat,
												task_args=[source]))

			handls.append(HandlerStat(inqu=qu,const=const))

			lev = parabam.Leviathon(procrs,handls,100000)
			lev.run()

			return out_paths

	def create_out_paths(self,outmode,output_group,user_structures):
		out_paths = {}
		for source in output_group:
			out_paths[source] = {}
			for name,struc in user_structures.items():
				if not struc.struc_type == np.ndarray:
					header,out_path = self.format_path(outmode,source,name)
					out_paths[source][name] = out_path
					with open(out_path,"w") as out_object:
						out_object.write(header)
				else:
					out_paths[source] = {name:""}
		return out_paths

	def format_path(self,mode,source,name):
		if mode == "a":
			header = "Sample,%s\n" % (self.name,)
			out_path = "./%s.csv" % (self.name,)
		elif mode == "s":
			header = "Analysis,Value\n"
			out_path = "./%s.csv" % (source,)

		return (header,out_path)

	def create_structures(self,user_structures_blueprint):
		#(data,store_method,)
		user_structures = {}
		class_to_type_map = {int:NumericStructure,float:NumericStructure,np.ndarray:ArrayStructure}
		for name,definition in user_structures_blueprint.items():
			definition["name"] = name
			definition["struc_type"] = type(definition["data"])
			user_structures[name] = class_to_type_map[definition["struc_type"]](**definition)
 		return user_structures

	def getParser(self):
		#argparse imported in ./interface/parabam 
		parser = ut.default_parser()
		parser.add_argument('--outmode', choices=['s','a'],default='s',
			help='Indicate whether data grouped by sample or analysis:\
			[s]ample: group output by sample,\
			[a]nalysis: group output by analysis')
		return parser 


#...happily ever after