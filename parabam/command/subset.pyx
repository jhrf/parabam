import pysam
import parabam
import time
import sys
import os
import gc
import shutil

from itertools import izip
from multiprocessing import Queue

class SubsetCore(object):

    def __init__(self,object constants,**kwargs):
        self.constants = constants

    def __pre_run_routine__(self,iterator,**kwargs):
        super(SubsetCore,self).__pre_run_routine__(iterator)
        self._user_subsets = self.constants.user_subsets
        self._counts = {}
        self._temp_paths = {}
        self._temp_objects = {}

        for subset in self._user_subsets:
            self._temp_paths[subset] = self.__get_temp_path__(subset)
            self._temp_objects[subset] = self.__get_temp_object__(self._temp_paths[subset])
            self._counts[subset] = 0

    def __get_results__(self):
        results = {}
        results["temp_paths"] = self._temp_paths
        results["counts"] = self._counts
        results["system"] = self._system

        return results

    def __write_to_subset_bam__(self,subset_type,read):
        self._counts[subset_type] += 1
        self._temp_objects[subset_type].write(read)

    def __post_run_routine__(self,**kwargs):
        super(SubsetCore,self).__post_run_routine__()
        for subset,fileobj in self._temp_objects.items():
            fileobj.close()


class Task(SubsetCore,parabam.command.Task):

    def __init__(self,parent_bam,inqu,outqu,task_size,constants,**kwargs):
        
        parabam.command.Task.__init__(self,parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    task_size=task_size,
                                    constants=constants)
        SubsetCore.__init__(self,constants)

    def __handle_engine_output__(self,engine_output,read):             
        if type(engine_output) == bool:
            if engine_output:
                self.__write_to_subset_bam__(self._constants.user_subsets[0],read)          
        elif type(engine_output) == list:
            for subset,cur_read in engine_output:
                self.__write_to_subset_bam__(subset,cur_read)
        elif type(engine_output) == None:
            pass

class PairTask(SubsetCore,parabam.command.PairTask):
    def __init__(self,parent_bam,inqu,outqu,task_size,constants):
        
        parabam.command.PairTask.__init__(self,parent_bam=parent_bam,
                                                inqu=inqu,
                                                outqu=outqu,
                                                task_size=task_size,
                                                constants=constants)
        SubsetCore.__init__(self,constants)

    def __handle_engine_output__(self,engine_output,read):
        for subset,cur_read in engine_output:
            self.__write_to_subset_bam__(subset,cur_read)

class Handler(parabam.command.Handler):

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object constants,object pause_qus,dict out_qu_dict):
        
        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                                     inqu=inqu,constants=constants,pause_qus=pause_qus,
                                     out_qu_dict=out_qu_dict)

        self._user_subsets = constants.user_subsets
        for subset in self._user_subsets:
            self._stage_stores[subset] = []

    def __new_package_action__(self,new_package):
        super(Handler,self).__new_package_action__(new_package)
        results = new_package.results
        for subset in self._user_subsets:
            self._stage_stores[subset].append((results["counts"][subset],
                                               results["temp_paths"][subset],))

    def __periodic_action__(self,iterations):
        super(Handler,self).__periodic_action__(iterations)

        for subset in self._user_subsets:
            if self.__test_stage_store__(subset):
                self.__add_merge_task__(results=self._stage_stores[subset],subset_type=subset)
                self.__update_stage_store__(subset)

    def __handler_exit__(self, **kwargs):
        super(Handler,self).__handler_exit__()
        #Kill the handlers handler

        for subset in self._user_subsets:
            #merge task 
            self.__add_merge_task__(self._stage_stores[subset],subset)
        self.__write_counts_csv__()
    
    def __add_merge_task__(self,results,subset_type):
        res = parabam.merger.MergePackage(results=results,subset_type=subset_type)
        self._out_qu_dict["merge"].put(res)

    def __write_counts_csv__(self):
        if self._constants.output_counts:
            count_path = "./subset_counts_csv_%d.csv" % (time.time(),)
            self.__standard_output__("[Status] Outputting counts to following file: %s" % (count_path,))

            with open(count_path, "w") as count_file:
                key_order,header = self.__get_count_header__()
                count_file.write(header)
                for line in self.__generate_count_line__(key_order):
                    count_file.write(line)

    def __generate_count_line__(self,key_order):
        for source,stats in self._stats.items():
            line = ""
            line += "%s" % (source,)
            for key in key_order:
                line += ",%s" % (stats[key],)
            line += "\n"
            yield line

    def __get_count_header__(self):
        header = "Sample"
        keys = []
        for source, stats in self._stats.items():
            for key, val in stats.items():
                header += ","
                header += key
                keys.append(key)
            break
        header += "\n"
        return keys, header

class Interface(parabam.command.Interface):
    """The interface to parabam subset.
    Users will primarily make use of the ``run`` function."""

    def __init__(self,temp_dir):
        super(Interface,self).__init__(temp_dir)

    def run_cmd(self,parser):
        cmd_args = parser.parse_args()

        module,user_engine,user_constants = self.__get_module_and_vitals__(cmd_args.instruc)

        if hasattr(module,"get_subset_types"):
            user_subsets = module.get_subset_types()
        else:
            user_subsets = ["subset"]

        self.run(
            input_bams=cmd_args.input,
            total_procs = cmd_args.p,
            task_size = cmd_args.s,
            verbose= cmd_args.v,
            user_subsets= user_subsets,
            user_constants = user_constants,
            user_engine = user_engine,
            fetch_region = cmd_args.region,
            reader_n = cmd_args.f,
            pair_process=cmd_args.pair,
            include_duplicates=cmd_args.d,
            debug = cmd_args.debug,
            ensure_unique_output=cmd_args.u,
            output_counts=cmd_args.counts
            )
    
    def run(self,input_bams,total_procs,task_size,user_constants,user_engine,
            user_subsets,reader_n = 2,fetch_region=None,
            keep_in_temp=False,verbose=0,pair_process=False,
            include_duplicates=True,debug=False,
            ensure_unique_output=False,output_counts=False):

        ''' Docstring! '''
        args = dict(locals())
        del args["self"]
        return super(Interface,self).run(**args)

    def __get_queue_names__(self,pair_process,**kwargs):
        queues = ["merge","main"]
        return queues

    def __get_destroy_handler_order__(self):
        return [Handler,parabam.merger.Handler]

    def __get_handler_bundle__(self,**kwargs):
        handler_bundle = {}

        #queues transformed by leviathon
        handler_bundle[Handler] = {"inqu":"main", 
                                   "out_qu_dict":["merge"]}

        handler_bundle[parabam.merger.Handler] = {"inqu":"merge",
                                                  "out_qu_dict":[]}
        return handler_bundle

    def __get_output_paths__(self,input_path,user_subsets,ensure_unique_output,**kwargs):
        output_paths = {input_path:{}}
        for salt,subset in enumerate(user_subsets):
            output_paths[input_path][subset] = self.__get_path__(input_path,
                                                                 subset,
                                                                 ensure_unique_output,
                                                                 salt)
        return output_paths

    def __get_path__(self,input_path,subset,ensure_unique_output,salt):
        head,ext = os.path.splitext(input_path)
        head,tail = os.path.split(head)
        unique = ""

        if ensure_unique_output:
            unique = "_%d%d" % (time.time(),salt)

        return os.path.join(".",self._temp_dir, "%s_%s%s%s" % (tail,subset,unique,ext))
            
    def __get_task_class__(self,pair_process,**kwargs):
        if pair_process:
            return PairTask
        else:
            return Task

    def __get_extra_const_args__(self,**kwargs):
        args = {}
        args["user_subsets"] = kwargs["user_subsets"]
        args["output_counts"] = kwargs["output_counts"]
        return args
        
    def get_parser(self):
        parser = self.default_parser()

        parser.add_argument('-u',action="store_true",default=False,
            help="The files will be appended with a unique identifier. In the format <input_name>_<unique_id>_<subset_type>.bam")
        parser.add_argument('--counts',action="store_true",default=False,
            help="The amount of reads in each subset will be output as a .CSV alongside the .BAM file.")

        return parser 

#...happily ever after
