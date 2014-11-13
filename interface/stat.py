import pysam
import parabam
import shutil
import time
import sys
import os
import copy
import numpy as np

from multiprocessing import Queue
from abc import ABCMeta, abstractmethod

class TaskStat(parabam.Task):

	__metaclass__ = ABCMeta

	def __init__(self,task_set,outqu,curproc,destroy,const,source):
		super(TaskStat, self).__init__(task_set=task_set,
										outqu=outqu,
										curproc=curproc,
										destroy=destroy,
										const=const)
		self._source = source
		self._master_file_path = const.master_file_path[self._source]
		self._counts = {}
		self._local_structures = {}

	def __generate_results__(self):
		master = pysam.Samfile(self._master_file_path,"rb")

		for name,struct in self.const.user_structures.items():
			self._local_structures[name] = struct.empty_clone()
			self._counts[name] = 0

		self.__handle_task_set__(self._task_set,master)

		master.close() #Close the master bamfile

		results = {}
		results["source"] = self._source
		results["structures"] = self.__unpack_structures__(self._local_structures)
		results["counts"] = self._counts

		return results

	def __unpack_structures__(self,structures):
		unpacked = []
		for name,struc in structures.items():
			unpacked.append( (name,struc.data) )

		return unpacked

	def __handle_task_set__(self,task_set,master):
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
		super(HandlerStat,self).__init__(inqu,const,destroy_limit=len(const.sources))

		self._sources = const.sources
		self._final_structures = {}

		for source in self._sources:
			self._final_structures[source] = {}
			for name,struc in self.const.user_structures.items():
				self._final_structures[source][struc.name] = struc.empty_clone() 

	def __new_package_action__(self,new_package,**kwargs):
		results = new_package.results
		source = results["source"]
		self.__auto_handle__(results,source)
		for name,data in results["structures"]:
			final_struc = self._final_structures[source][name]
			final_struc.merge(data)

	def __periodic_action__(self,iterations):
		pass

	def __handler_exit__(self,**kwargs):
		const = self.const
		if const.outmode == "d":
			for source in self._sources:
				source_structures = self._final_structures[source]
				data_str = self.__get_data_str_from_names__(const.analysis_names,source_structures)
				with open(const.output_paths,"a") as out_object:
					out_object.write("%s%s\n" % (source,data_str))

				#Arrays can't be squished into unified output, so create a unique file path	
				for name,struc in source_structures.items():
					if struc.struc_type == np.ndarray:
						array_path = const.output_paths.replace("%s/" % (const.temp_dir,),"./%s_" % (struc.name,))
						struc.write_to_csv(array_path,source,const.outmode)
		else:
			for source in self._sources:
				for name,struc in self._final_structures[source].items():
					struc.write_to_csv(const.output_paths,source,const.outmode)

	def __get_data_str_from_names__(self,names,user_structures):
		data_str = ""
		for name in names:
			cur_data = user_structures[name].data
			data_str += ",%.3f" % (cur_data,)
		return data_str

class ProcessorStat(parabam.Processor):

	def __init__(self,outqu,const,TaskClass,task_args):
		super(ProcessorStat,self).__init__(outqu,const,TaskClass,task_args)
		self._source = task_args[0] #Defined in the run function within Interface

	def __get_master_bam__(self,master_file_path):
		return pysam.Samfile(master_file_path[self._source],"rb")

	def __add_to_collection__(self,master,alig,collection):
		collection.append(alig)

	def __pre_processor__(self,master_file_path):
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
		super(ArrayStructure,self).__init__(name,struc_type,store_method,data)

	def empty_clone(self):
		return ArrayStructure(self.name,self.struc_type,self.store_method,self.org_data)

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

	def write_to_csv(self,out_path,source,mode):
		np.savetxt(out_path,self.data,fmt="%.3f",delimiter=",")

class Interface(parabam.Interface):

	def __init__(self,temp_dir,exe_dir):
		super(Interface,self).__init__(temp_dir,exe_dir)
	
	def run_cmd(self,parser):

		cmd_args = parser.parse_args()

		verbose = cmd_args.v
		module = __import__(cmd_args.instruc, fromlist=[''])
		user_engine = module.engine
		user_constants = {}
		user_structures_blueprint = {}
		module.set_constants(user_constants)
		module.set_structures(user_structures_blueprint)

		self.run(
			input_bams=cmd_args.input,
			output_path = cmd_args.output,
			proc= cmd_args.p,
			chunk= cmd_args.c,
			verbose= verbose,
			user_constants = user_constants,
			user_engine = user_engine,
			user_struc_blueprint = user_structures_blueprint,
			engine_is_class = False,
			outmode = cmd_args.outmode)
	
	def run(self,input_bams,output_path,proc,chunk,user_constants,user_engine,
			user_struc_blueprint,outmode,multi=4,verbose=False,engine_is_class=False):

		user_structures = self.__create_structures__(user_struc_blueprint)
		super_sources = [ self.__get_basename__(b) for b in input_bams ]
		if not output_path:
			output_path = "%s/parabam_stat_%d.csv" % (self._temp_dir,int(time.time()),)

		analysis_names = self.__get_non_array_names__(user_structures)
		if outmode == "d":
			header = "Sample,%s\n" % (",".join(analysis_names),)
			with open(output_path,"w") as out_obj:
				out_obj.write(header)
		master_out_paths = {}
		for input_group,output_group in self.__get_group__(input_bams,super_sources,multi=multi):
			
			task_qu = Queue()

			master_file_path = {}

			for master,source in zip(input_group,output_group):
				master_file_path[source] = master
					
			procrs = []
			handls = []

			if outmode == "d":
				output_paths = output_path
				master_output_paths = output_path
			else:
				output_paths = self.__create_output_paths__(outmode,output_group,user_structures)
				master_output_paths.update(output_paths)

			const = parabam.Const(output_paths=output_paths,temp_dir=self._temp_dir,
								master_file_path=master_file_path,
								chunk=chunk,proc=(proc // len(input_group)),
								verbose=verbose,thresh=0,
								sources=output_group,
								exe_dir=self._exe_dir,
								user_constants=user_constants,
								user_structures=user_structures,
								user_engine=user_engine,
								outmode=outmode,
								analysis_names = analysis_names)

			for source in output_group:
				if engine_is_class:
					if not issubclass(user_engine,TaskStat):
						raise Exception("[ERROR]\tThe class provided to parabam must be a subclass of\n"\
										"\tparabam.interface.subset.TaskStat. Please consult the parabam manual.")
					cur_args = [source]
					procrs.append(ProcessorStat(outqu=task_qu,
											const=const,
											TaskClass=user_engine,
											task_args=cur_args))
				else:
					procrs.append(ProcessorStat(outqu=task_qu,
												const=const,
												TaskClass=TaskStat,
												task_args=[source]))

			handls.append(HandlerStat(inqu=task_qu,const=const))

			lev = parabam.Leviathon(procrs,handls,100000)
			lev.run()

		shutil.move(output_path,"./")
		return master_output_paths #Depending on outmode this output a dict or string

	def __create_output_paths__(self,outmode,output_group,user_structures):
		out_paths = {}
		for source in output_group:
			out_paths[source] = {}
			for name,struc in user_structures.items():
				if not struc.struc_type == np.ndarray:
					header,out_path = self.__format_path__(outmode,source,name)
					out_paths[source][name] = out_path
					with open(out_path,"w") as out_object:
						out_object.write(header)
				else:
					out_paths[source] = {name:""}
		return out_paths

	def __format_path__(self,mode,source,name):
		if mode == "a":
			header = "Sample,%s\n" % (self.name,)
			out_path = "./%s.csv" % (self.name,)
		elif mode == "s":
			header = "Analysis,Value\n"
			out_path = "./%s.csv" % (source,)
		return (header,out_path)

	def __get_non_array_names__(self,user_structures):
		names = []
		for name,struc in user_structures.items():
			if not struc.struc_type == np.ndarray:
				names.append(name)
		names.sort()
		return names

	def __create_structures__(self,user_structures_blueprint):
		#(data,store_method,)
		user_structures = {}
		class_to_type_map = {int:NumericStructure,float:NumericStructure,np.ndarray:ArrayStructure}
		for name,definition in user_structures_blueprint.items():
			definition["name"] = name
			definition["struc_type"] = type(definition["data"])
			user_structures[name] = class_to_type_map[definition["struc_type"]](**definition)
 		return user_structures

	def get_parser(self):
		#argparse imported in ./interface/parabam 
		parser = self.default_parser()

		parser.add_argument('--instruc','-i',metavar='INSTRUCTION',required=True
			,help='The instruction file, written in python, that we wish'\
			'to carry out on the input BAM.')
		parser.add_argument('--input','-b',metavar='INPUT', nargs='+',required=True
			,help='The file(s) we wish to operate on. Multipe entries should be separated by a single space')
		parser.add_argument('--output','-o',metavar='OUTPUT', nargs='?',required=False
		,help='Specify a name for the output CSV file. Only used with default `outmode`.\
			If this argument is not supplied, the output will take the following form:\
			parabam_stat_[UNIX_TIME].csv')
		parser.add_argument('--outmode', choices=['d','s','a'],default='d',
			help='Indicate whether data grouped by sample or analysis:\
			[d]efault: a csv with a column for each analysis and row for each sample\
			[s]ample: csv for each sample,\
			[a]nalysis: csv for each analysis')
		return parser 


#...happily ever after
