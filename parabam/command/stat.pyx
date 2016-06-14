#Once upon a time...
import pysam
import parabam
import shutil
import time
import sys
import os
import copy
import numpy as np

from itertools import izip
from multiprocessing import Queue
from abc import ABCMeta, abstractmethod

class StatCore(object):

    def __init__(self,object constants):
        self._constants
        self._counts = {}
        self._local_structures = {}
        self._system = {}

    def __pre_run_routine__(self,iterator,**kwargs):
        super(StatCore,self).__pre_run_routine__(iterator)
        for name,structure in self._constants.user_structures.items():
            self._local_structures[name] = structure.empty_clone()
            self._counts[name] = 0

    def __get_results__(self):
        results = {}
        results["structures"] = \
                    self.__unpack_structures__(self._local_structures)
        results["counts"] = self._counts
        results["system"] = self._system
        return results

    def __post_run_routine__(self,**kwargs):
        super(StatCore,self).__post_run_routine__()
        pass

    def __handle_rule_output__(self,rule_output,read):
        if rule_output:#Allows return False
            local_structures = self._local_structures
            for name,result in rule_output.items():
                self._counts[name] += 1
                local_structures[name].add(result)

    def __unpack_structures__(self,structures):
        unpacked = []
        for name,struc in structures.items():
            unpacked.append( (name,struc.data) )
        return unpacked

class Task(StatCore,parabam.command.Task):

    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants):
        
        parabam.command.Task.__init__(self,parent_bam=parent_bam,
                                            inqu=inqu,
                                            outqu=outqu,
                                            statusqu=statusqu,
                                            task_size=task_size,
                                            constants=constants)
        StatCore.__init__(self,constants)
        
class PairTask(StatCore,parabam.command.PairTask):

    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants):
        
        parabam.command.PairTask.__init__(self,parent_bam=parent_bam,
                                            inqu=inqu,
                                            outqu=outqu,
                                            statusqu=statusqu,
                                            task_size=task_size,
                                            constants=constants)
        StatCore.__init__(self,constants)


class Handler(parabam.command.Handler):

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object constants,object pause_qus,dict out_qu_dict,object report=True):

        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                                     inqu=inqu,constants=constants,pause_qus=pause_qus,
                                     out_qu_dict=out_qu_dict)

        self._final_structures = {}
        for name,struc in self._constants.user_structures.items():
            self._final_structures[struc.name] = struc.empty_clone()

    def __new_package_action__(self,new_package,**kwargs):
        super(Handler,self).__new_package_action__(new_package)
        results = new_package.results
        for name,data in results["structures"]:
            final_struc = self._final_structures[name]
            final_struc.merge(data)

    def __handler_exit__(self):
        super(Handler,self).__handler_exit__()
        constants = self._constants

        #Append to global csv
        if constants.numeric_names: #Check that there are global analyses
            data_str = \
                self.__get_data_str_from_names__(constants.numeric_names,
                                                  self._final_structures)
            
            print self._output_paths
            with open(self._output_paths["global"]["stats"],"a") as out_object:
                out_object.write("%s%s\n" % \
                    (self._parent_bam.filename,data_str))
        
        #Output non global data
        for name,structure in self._final_structures.items():
            if structure.struc_type == np.ndarray or \
                structure.struc_type == dict:
                structure.write_to_csv(\
                    self._output_paths[self._parent_bam.filename][name])

    def __get_data_str_from_names__(self,names,user_structures):
        data_str = ""
        for name in names:
            if user_structures[name].struc_type == np.ndarray:
                continue
            cur_data = user_structures[name].data
            data_str += ",%.5f" % (cur_data,)
        return data_str

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

class NumericStructure(UserStructure):
    def __init__(self,name,struc_type,store_method,data,log_scaling=False):
        super(NumericStructure,self).__init__(name,struc_type,store_method,data)
        self.log_scaling = log_scaling

        if store_method == "min":
            self.data = float('inf')
            self.org_data = copy.copy(float('inf'))

    def empty_clone(self):
        return NumericStructure(self.name,
                                self.struc_type,
                                self.store_method,
                                self.org_data)

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

#TODO: This mode doesn't work at all. Probably something to do with 
#      creating an empty clone. Counts are inflated.
class CounterStructure(UserStructure):
    def __init__(self,name,struc_type,store_method,data):
        super(CounterStructure,self).__init__(name,struc_type,store_method,data)

    def empty_clone(self):
        return CounterStructure(self.name,self.struc_type,self.store_method,self.org_data)

    def add_cumu(self,result):
        for key,value in result.items():
            try:
                self.data[key] += value
            except KeyError:
                self.data[key] = value

    def add_max(self,result):
        for key,value in result.items():
            try:
                self.data[key] = max([self.data[key],value])
            except KeyError:
                self.data[key] = value

    def add_min(self,result):
        for key,value in result.items():
            try:
                self.data[key] = min([self.data[key],value])
            except KeyError:
                self.data[key] = value

    def merge_cumu(self,result):
        self.add_cumu(result)
        del result

    def merge_max(self,result):
        self.add_max(result)
        del result

    def merge_min(self,result):
        self.add_min(result)
        del result

    def write_to_csv(self,out_path):

        with open(out_path,"w") as out_file:
            for key,value in self.data.items():
                out_str = "%s,%.5f\n" % (key,value,)
                out_file.write(out_str)

class ArrayStructure(UserStructure):
    def __init__(self,name,struc_type,store_method,data):
        super(ArrayStructure,self).__init__(name,struc_type,store_method,data)

        if self.store_method == "vstack":
            self.seen = 0
            self.add = self.add_vstack
            self.merge = self.merge_vstack

    def empty_clone(self):
        return ArrayStructure(self.name,self.struc_type,self.store_method,self.org_data)

    def add_max(self,result,coords):
        existing = self.data[coords]
        self.data[coords] = self.max_decision(result,existing)

    def add_min(self,result,coords):
        existing = self.data[coords]
        self.data[coords] = self.min_decision(result,existing)

    def add_cumu(self,result):
        self.data = np.add(self.data,result)

    def add_vstack(self,result):
        if self.seen == 0:
            self.data = result
            self.seen += 1
        else:
            self.data = np.vstack((self.data,result))

    def merge_max(self,result):
        self.data = np.maximum(self.data,result)
        del result

    def merge_min(self,result):
        self.data = np.minimum(self.data,result)
        del result

    def merge_cumu(self,result):
        self.data = np.add(self.data,result)
        del result

    def merge_vstack(self,result):
        self.add_vstack(result)
        del result

    def write_to_csv(self,out_path):

        format = []
        for x in self.data[0,:]:
            type_of_x = type(x) 
            if type_of_x == str or type_of_x == np.string_:
                format.append("%s")
            else:
                format.append("%.5f")

        np.savetxt(out_path,self.data,fmt=",".join(format),delimiter=",")

class Stat(parabam.command.Interface):

    def __init__(self,**kwargs):
        super(Stat,self).__init__(instance_name = "parabam stat", **kwargs)
    
    def run_cmd(self):
        module,user_rule,user_constants = \
                self.__get_module_and_vitals__(self.cmd_args.rule)

        user_struc_blueprint = {}
        module.get_blueprints(user_struc_blueprint)

        self.run(input_paths=self.cmd_args.input,
                  user_constants = user_constants,
                  user_rule = user_rule,
                  user_struc_blueprint = user_struc_blueprint,
                  fetch_region = self.cmd_args.region)

    def run(self,input_paths,
                  user_constants,
                  user_rule,
                  user_struc_blueprint,
                  user_specified_outpath=None,
                  fetch_region=None,
                  **kwargs):
                  
        ''' Docstring! '''
        args = dict(locals())
        del args["self"]

        #Prepare state structures and insert to args
        #kwargs are later used to construct the Constant file
        #passed to all the fileprocessors and handlers
        user_structures = self.__create_structures__(user_struc_blueprint)
        numeric_names = self.__get_numeric_names__(user_structures)

        args["user_structures"] = user_structures
        args["numeric_names"] = numeric_names

        results = super(Stat,self).run(**args)
        return results

    def __get_destroy_handler_order__(self):
        return [Handler]

    def __get_queue_names__(self,**kwargs):
        return ["main"]

    def __get_handler_bundle__(self,**kwargs):
        handler_bundle = { Handler: {"inqu":"main","out_qu_dict":[]}}
        return handler_bundle

    def __instalise_final_output__(self,numeric_names,
                                        user_specified_outpath,
                                        **kwargs):
        final_output = {}

        if len(numeric_names) > 0:
            global_filename = \
                        self.__get_global_output_path__(user_specified_outpath)            
            self.__create_global_output_file__(global_filename,numeric_names)

            final_output["global"] = {"stats": global_filename}

        return final_output

    def __get_global_output_path__(self,user_specified_outpath):
        
        if user_specified_outpath is None:
            global_filename = \
                    os.path.join(self.temp_dir,"parabam_stat_%d_%d.csv"\
                                    % (time.time(),os.getpid()))
        else:
            global_filename = user_specified_outpath

        return global_filename

    def __get_output_paths__(self,input_path,
                             final_output_paths,
                             user_structures,
                             **kwargs):
        output_paths = {input_path:{}}

        for name,structure in user_structures.items():
            if not issubclass(structure.__class__,NumericStructure):
                path_id,ext = os.path.splitext(os.path.basename(input_path))
                csv_path = "%s_%s.csv" % (path_id,name,)
                output_paths[input_path][name] =\
                            os.path.join(".",self.temp_dir,csv_path)

        if "global" in final_output_paths.keys():
            output_paths["global"] = final_output_paths["global"] 
        return output_paths

    def __get_task_class__(self,**kwargs):
        if self.pair_process:
            return PairTask
        else:
            return Task

    def __get_queues__(self,object constants,**kwargs):
        queues = {"main":Queue()}
        return queues

    def __create_global_output_file__(self,output_path,numeric_names):
        header = "Sample,%s\n" % (",".join(numeric_names),)
        with open(output_path,"w") as out_obj:
            out_obj.write(header)

    def __get_numeric_names__(self,user_structures):
        numeric_analysis = []
        for name,structure in user_structures.items():
            if issubclass(structure.__class__,NumericStructure):
                numeric_analysis.append(name)
        return sorted(numeric_analysis)

    def __create_structures__(self,user_struc_blueprint):
        user_structures = {}
        class_to_type_map = {int:NumericStructure,
                             float:NumericStructure,
                             np.ndarray:ArrayStructure,
                             dict:CounterStructure}

        for name,definition in user_struc_blueprint.items():
            definition["name"] = name
            definition["struc_type"] = type(definition["data"])
            user_structures[name] = \
                    class_to_type_map[definition["struc_type"]](**definition)
        return user_structures

    def get_parser(self):
        #argparse imported in ./interface/parabam 
        parser = self.default_parser()

        parser.add_argument('--output','-o',
                                metavar='OUTPUT', 
                                nargs='?',
                                required=False
        ,help="Specify a name for the output CSV file. If this argument is \n"\
            "not supplied, the output will take the following form:\n"\
            "parabam_stat_[UNIX_TIME].csv")

        return parser

#...happily ever after
