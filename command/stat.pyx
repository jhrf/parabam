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

    def __init__(self,object const,source):
        self._source = source
        self._counts = {}
        self._local_structures = {}
        self._system = {}

    def __generate_results__(self):
        for name,structure in self.const.user_structures.items():
            self._local_structures[name] = structure.empty_clone()
            self._counts[name] = 0

        self.__process_task_set__(self._task_set)

        results = {}
        results["structures"] = self.__unpack_structures__(self._local_structures)
        results["counts"] = self._counts
        results["system"] = self._system

        return results

    def __handle_engine_output__(self,engine_output,read):
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

    def __init__(self, object task_set, object outqu, object curproc,object destroy,
                 object parent_bam, object const, str parent_class,str source):
        
        StatCore.__init__(self,const,source)
        parabam.command.Task.__init__(self,task_set=task_set,
                                        outqu=outqu,
                                        curproc=curproc,
                                        destroy=destroy,
                                        const=const,
                                        parent_class = parent_class,
                                        parent_bam = parent_bam)
        
class PairTask(StatCore,parabam.command.PairTask):

    def __init__(self, object task_set, object outqu, object curproc,object destroy,
                 object parent_bam,object const, str parent_class,str source):
        
        StatCore.__init__(self,const,source)
        parabam.command.PairTask.__init__(self,task_set=task_set,
                                        outqu=outqu,
                                        curproc=curproc,
                                        destroy=destroy,
                                        const=const,
                                        parent_class = parent_class,
                                        parent_bam = parent_bam)


class Handler(parabam.command.Handler):

    def __init__(self,object inqu,object const,int destroy_limit,dict out_qu_dict):
        super(Handler,self).__init__(inqu,const,destroy_limit=destroy_limit,out_qu_dict=out_qu_dict)
        self._final_structures = {}

        for source in self._sources:
            self._final_structures[source] = {}
            for name,struc in self.const.user_structures.items():
                self._final_structures[source][struc.name] = struc.empty_clone()

    def __new_package_action__(self,new_package,**kwargs):
        super(Handler,self).__new_package_action__(new_package)
        results = new_package.results
        source = results["source"]
        for name,data in results["structures"]:
            final_struc = self._final_structures[source][name]
            final_struc.merge(data)

    def __handler_exit__(self):
        super(Handler,self).__handler_exit__()
        const = self.const

        for source in self._sources:
            source_structures = self._final_structures[source]
            data_str = self.__get_data_str_from_names__(const.analysis_names,source_structures)

            with open(const.output_paths["global"][0],"a") as out_object:
                out_object.write("%s%s\n" % (source,data_str))

            #Arrays can't be squished into unified output, so create a unique file path    
            for name,struc in source_structures.items():
                if struc.struc_type == np.ndarray:
                    array_path = "%s_%s.csv" % (source,struc.name)
                    struc.write_to_csv(array_path,source,const.outmode)

    def __get_data_str_from_names__(self,names,user_structures):
        data_str = ""
        for name in names:
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
            self.data = float('inf')
            self.org_data = copy.copy(float('inf'))
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
        existing = self.data[coords]
        self.data[coords] = self.max_decision(result,existing)

    def add_min(self,result,coords):
        existing = self.data[coords]
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

class Interface(parabam.command.Interface):

    def __init__(self,temp_dir):
        super(Interface,self).__init__(temp_dir)
    
    def run_cmd(self,parser):

        cmd_args = parser.parse_args()

        module,user_engine,user_constants = self.__get_module_and_vitals__(cmd_args.instruc)
        user_struc_blueprint = {}
        module.set_structures(user_struc_blueprint)

        self.run(
            input_bams=cmd_args.input,
            user_specified_outpath=cmd_args.output,
            proc=cmd_args.p,
            chunk=cmd_args.c,
            verbose= cmd_args.v,
            user_constants = user_constants,
            user_engine = user_engine,
            user_struc_blueprint = user_struc_blueprint,
            fetch_region=cmd_args.region,
            side_by_side=cmd_args.s,
            pair_process=cmd_args.pair)
    
    def run(self,input_bams,proc,chunk,user_constants,user_engine,
            user_struc_blueprint,user_specified_outpath=None,
            fetch_region=None,side_by_side=2,keep_in_temp=False,
            engine_is_class=False,verbose=0,pair_process=False,
            include_duplicates=True,ensure_unique_output=False,
            debug=False,input_is_sam=False):

        ''' Docstring! '''
        args = dict(locals())
        del args["self"]
        return super(Interface,self).run(**args)

    def __get_processor_bundle__(self,queues,object const,task_class,pair_process,debug,**kwargs):
        processor_bundle = {}

        processor_bundle["class"] = parabam.command.Processor

        for i,source in enumerate(const.sources):
            processor_bundle[source] = {"source":source,
                                        "outqu":queues["main"],
                                        "const":const,
                                        "TaskClass":task_class,
                                        "task_args":[source],
                                        "debug":debug}
        return processor_bundle


    def __get_handler_bundle__(self,queues,object const,task_class,**kwargs):
        handler_bundle = {Handler: {"inqu":queues["main"],
                                    "const":const,
                                    "destroy_limit":len(const.sources),
                                    "out_qu_dict":{}} }
        return handler_bundle

    def __get_master_file_path__(self,sources,input_paths,**kwargs):
        master_file_path = {}
        for path,source in izip(input_paths,sources):
                master_file_path[source] = path
        return master_file_path

    def __get_output_paths__(self,sources,input_paths,user_specified_outpath,
                             user_struc_blueprint,**kwargs):
        output_paths = {}

        if not user_specified_outpath:
            output_paths["global"] = [os.path.join(".",self._temp_dir,
                                                   "parabam_stat_%d.csv" % (time.time(),))]
        else:
            output_paths["global"] = [user_specified_outpath]
        analysis_names = self.__get_non_array_names__(user_struc_blueprint)
        self.__create_output_files__(output_paths["global"][0],analysis_names)

        for name,blueprint in user_struc_blueprint.items():
            if blueprint["data"] == np.ndarray:
                try:
                    path_list = output_paths[source]
                except KeyError:
                    output_paths[source] = []
                    path_list = output_paths[source]

                for source in sources:
                    filename = "%s_%s_array.csv" % (source,name,)
                    path_list.append(os.path.join(".",self._temp_dir,filename))

        return output_paths

    def __get_const_args__(self,**kwargs):
        args = super(Interface,self).__get_const_args__(**kwargs)
        user_structures = self.__create_structures__(kwargs["user_struc_blueprint"])
        analysis_names = self.__get_non_array_names__(kwargs["user_struc_blueprint"])
        args["user_structures"] = user_structures
        args["analysis_names"] = analysis_names
        return args

    def __get_task_class__(self,pair_process,user_engine,engine_is_class,**kwargs):
        if pair_process:
            base_class = PairTask
        else:
            base_class = Task

        if engine_is_class:
            if not pair_process and not issubclass(user_engine,Task):
                raise_exception = True
            elif pair_process and not issubclass(user_engine,PairTask):
                raise_exception = True
            else:
                raise_exception = False
                
            if raise_exception:
                #TODO: Untested exception here. Class name might be ABCMeta also
                raise Exception("[ERROR] User engine class must inherit %s\n" \
                    % (base_class.__class__,))
                sys.exit(1)
            else:
                task_class = user_engine
        else:
            task_class = base_class

        return task_class

    def __get_queues__(self,object const,**kwargs):
        queues = {"main":Queue()}
        return queues

    def __create_output_files__(self,output_path,analysis_names):
        header = "Sample,%s\n" % (",".join(analysis_names),)
        with open(output_path,"w") as out_obj:
            out_obj.write(header)

    def __get_non_array_names__(self,user_struc_blueprint):
        names = []
        for name,blueprint in user_struc_blueprint.items():
            if not blueprint["data"] == np.ndarray:
                names.append(name)
        names.sort()
        return names

    def __create_structures__(self,user_struc_blueprint):
        #(data,store_method,)
        user_structures = {}
        class_to_type_map = {int:NumericStructure,float:NumericStructure,np.ndarray:ArrayStructure}
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
