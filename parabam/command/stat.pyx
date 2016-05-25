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
        results["structures"] = self.__unpack_structures__(self._local_structures)
        results["counts"] = self._counts
        results["system"] = self._system
        return results

    def __post_run_routine__(self,**kwargs):
        super(StatCore,self).__post_run_routine__()
        pass

    def __handle_engine_output__(self,engine_output,read):
        if engine_output:#Allows return False
            local_structures = self._local_structures
            for name,package in engine_output.items():
                self._counts[name] += 1
                local_structures[name].add(**package)

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
        if constants.analysis_names: #Check that there are global analyses
            data_str = self.__get_data_str_from_names__(constants.analysis_names
                                                        ,self._final_structures)
            
            with open(self._output_paths["global"][0],"a") as out_object:
                out_object.write("%s%s\n" % \
                    (self._parent_bam.filename,data_str))
        
        #Output numpy arrays
        for name,structure in self._final_structures.items():
            if structure.struc_type == np.ndarray or \
                structure.struc_type == dict:
                structure.write_to_csv(self._output_paths[self._parent_bam.filename][name])

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

class Interface(parabam.command.Interface):

    def __init__(self,temp_dir):
        super(Interface,self).__init__(temp_dir)
    
    def run_cmd(self,parser):

        cmd_args = parser.parse_args()

        module,user_engine,user_constants = self.__get_module_and_vitals__(cmd_args.instruc)
        user_struc_blueprint = {}
        module.set_structures(user_struc_blueprint)

        self.run(input_bams=cmd_args.input,
            total_procs = cmd_args.p,
            task_size = cmd_args.s,
            verbose= cmd_args.v,
            reader_n = cmd_args.f,
            user_constants = user_constants,
            user_struc_blueprint = user_struc_blueprint,
            user_engine = user_engine,
            fetch_region = cmd_args.region,
            pair_process=cmd_args.pair,
            coord_process=cmd_args.coord,
            include_duplicates=cmd_args.d,
            debug = cmd_args.debug,
            announce = True)


    def run(self,input_bams,total_procs,task_size,user_constants,user_engine,
            user_struc_blueprint,user_specified_outpath=None,
            reader_n = 2,fetch_region=None,side_by_side=2,
            keep_in_temp=False,verbose=0,coord_process=False,
            pair_process=False,include_duplicates=True,
            debug=False,announce=False):

        ''' Docstring! '''
        args = dict(locals())
        del args["self"]

        if not verbose:
            announce = False
        self.__introduce__("parabam stat",announce)

        results = super(Interface,self).run(**args)

        self.__goodbye__("parabam stat",announce)
        return results

    def __get_destroy_handler_order__(self):
        return [Handler]

    def __get_queue_names__(self,**kwargs):
        return ["main"]

    def __get_extra_const_args__(self,user_struc_blueprint,**kwargs):
        user_structures = self.__create_structures__(user_struc_blueprint)
        analysis_names = self.__get_non_array_names__(user_struc_blueprint)

        return {"user_structures":user_structures,
                "analysis_names":analysis_names}

    def __get_handler_bundle__(self,**kwargs):
        handler_bundle = { Handler: {"inqu":"main",
                                    "out_qu_dict":[]}}
        return handler_bundle


    def __get_global_output_path__(self,user_struc_blueprint,user_specified_outpath,**kwargs):
        analysis_names = self.__get_non_array_names__(user_struc_blueprint)
        if analysis_names:
            if not user_specified_outpath:
                global_filename = os.path.join(".",self._temp_dir,"parabam_stat_%d_%d.csv"\
                                         % (time.time(),os.getpid()))
            else:
                global_filename = user_specified_outpath
            self.__create_output_files__(global_filename,analysis_names)
            return {"global": [global_filename]}
        else:
            return {"global": []}

    def __get_output_paths__(self,input_path,user_specified_outpath,
                             user_struc_blueprint,
                             **kwargs):
        output_paths = {input_path:{}}

        # TODO: See note on line 428. Related problem here
        for name,blueprint in user_struc_blueprint.items():
            if type(blueprint["data"]) == np.ndarray or \
                     type(blueprint["data"]) == dict:
                path_id,ext = os.path.splitext(os.path.basename(input_path))
                csv_path = "%s_%s.csv" % (path_id,name,)
                output_paths[input_path][name] = os.path.join(".",self._temp_dir,csv_path)
        return output_paths

    def __get_task_class__(self,pair_process,**kwargs):
        if pair_process:
            return PairTask
        else:
            return Task

    def __get_queues__(self,object constants,**kwargs):
        queues = {"main":Queue()}
        return queues

    def __create_output_files__(self,output_path,analysis_names):
        header = "Sample,%s\n" % (",".join(analysis_names),)
        with open(output_path,"w") as out_obj:
            out_obj.write(header)

    def __get_non_array_names__(self,user_struc_blueprint):
        names = []

        # TODO: dict type has been squeezed in here and as a result
        #       this function name no longer makes sense
        #
        #       consider refactoring so that numeric types are the
        #       exception rather than having not(numeric_types) as
        #       the exceptional case
        for name,blueprint in user_struc_blueprint.items():
            if not type(blueprint["data"]) == np.ndarray and \
                not type(blueprint["data"]) == dict:
                names.append(name)
        names.sort()
        return names

    def __create_structures__(self,user_struc_blueprint):
        user_structures = {}
        class_to_type_map = {int:NumericStructure,
                             float:NumericStructure,
                             np.ndarray:ArrayStructure,
                             dict:CounterStructure}

        for name,definition in user_struc_blueprint.items():
            definition["name"] = name
            definition["struc_type"] = type(definition["data"])
            user_structures[name] = class_to_type_map[definition["struc_type"]](**definition)
        return user_structures

    def get_parser(self):
        #argparse imported in ./interface/parabam 
        parser = self.default_parser()

        parser.add_argument('--output','-o',metavar='OUTPUT', nargs='?',required=False
        ,help="Specify a name for the output CSV file. Only used with default `outmode`.\n"\
            "If this argument is not supplied, the output will take the following form:\n"\
            "parabam_stat_[UNIX_TIME].csv")

        return parser

#...happily ever after
