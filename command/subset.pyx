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

    def __init__(self,object const):
        self._user_subsets = const.user_subsets
        self._counts = {}
        self._system = {}
        self._temp_paths = {}
        self._temp_objects = {}

    def __generate_results__(self):
        cdef dict temp_paths = self._temp_paths
        cdef dict temp_objects = self._temp_objects
        cdef dict counts = self._counts
        cdef dict system = self._system

        for subset in self._user_subsets:
            ext = self.__get_extension__(self.output_paths[subset])
            temp_paths[subset] = self.__get_temp_path__(subset,ext)
            temp_objects[subset] = self.__get_temp_object__(temp_paths[subset],ext)
            counts[subset] = 0

        self.__process_task_set__(self._task_set)

        for subset,file_object in temp_objects.items():
            file_object.close() #close all the other bams
        del self._temp_objects

        results = {}
        results["temp_paths"] = temp_paths
        results["counts"] = counts
        results["system"] = system

        return results

    def __write_to_subset_bam__(self,subset_type,read):
        self._counts[subset_type] += 1
        self._temp_objects[subset_type].write(read)

class Task(SubsetCore,parabam.command.Task):

    def __init__(self, object task_set, object outqu, object curproc,
                 object parent_bam, object const, str parent_class):
        SubsetCore.__init__(self,const)
        parabam.command.Task.__init__(self,task_set=task_set,
                                      outqu=outqu,
                                      curproc=curproc,
                                      parent_bam = parent_bam,
                                      const=const,
                                      parent_class=parent_class)

    def __handle_engine_output__(self,engine_output,read):
        subset_write = self.__write_to_subset_bam__
             
        if type(engine_output) == bool:
            if engine_output:
                subset_write(self.const.user_subsets[0],read)          
        elif type(engine_output) == list:
            for subset,cur_read in engine_output:
                subset_write(subset,cur_read)
        elif type(engine_output) == int:
            if not engine_output == -1:
                subset_write(self.const.user_subsets[engine_output],read)
        elif type(engine_output) == None:
            pass

class PairTask(SubsetCore,parabam.command.PairTask):
    def __init__(self, object task_set, object outqu, object curproc,
                 object parent_bam,object const,
                 str parent_class):
        
        SubsetCore.__init__(self,const)
        parabam.command.PairTask.__init__(self,task_set=task_set,
                        outqu=outqu,
                        curproc=curproc,
                        parent_bam = parent_bam,
                        const=const,
                        parent_class=parent_class)

    def __handle_engine_output__(self,engine_output,read):
        for subset,cur_read in engine_output:
            self.__write_to_subset_bam__(subset,cur_read)

class Handler(parabam.command.Handler):

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object const,object pause_qu,dict out_qu_dict):
        
        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                                     inqu=inqu,const=const,pause_qu=pause_qu,
                                     out_qu_dict=out_qu_dict)

        self._user_subsets = const.user_subsets
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
            self.__add_merge_task__(subset)
        self.__write_counts_csv__()
    
    def __add_merge_task__(self,results,subset_type):
        res = parabam.merger.MergePackage(results=results,subset_type=subset_type)
        self._out_qu_dict["merge"].put(res)

    def __write_counts_csv__(self):
        if self.const.output_counts:
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
            proc= cmd_args.p,
            chunk= cmd_args.c,
            verbose= cmd_args.v,
            user_subsets= user_subsets,
            user_constants = user_constants,
            user_engine = user_engine,
            engine_is_class = False,
            fetch_region = cmd_args.region,
            pair_process=cmd_args.pair,
            include_duplicates=cmd_args.d,
            side_by_side = cmd_args.s,
            debug = cmd_args.debug,
            ensure_unique_output=cmd_args.u,
            output_counts=cmd_args.counts
            )
    
    def run(self,input_bams,proc,chunk,user_constants,user_engine,
            user_subsets,fetch_region=None,side_by_side=2,
            keep_in_temp=False,engine_is_class=False,verbose=0,
            pair_process=False,include_duplicates=True,debug=False,
            ensure_unique_output=False,output_counts=False,
            input_is_sam=False):

        ''' Docstring! '''
        args = dict(locals())
        del args["self"]
        return super(Interface,self).run(**args)

    def __get_queues__(self,object const,pair_process,**kwargs):
        queues = ["merge","main"]
        if pair_process:
            queues = ["chaser"]

        return queues

    def __get_processor_bundle__(self,object const,task_class,debug,**kwargs): 
        processor_bundle = {"outqu":"main",
                            "const":const,
                            "TaskClass":task_class,
                            "task_args":{},
                            "debug":debug,
                            "class":parabam.core.Processor}

        return processor_bundle

    def __get_handler_bundle__(self,object const,task_class,pair_process,**kwargs):
        handler_bundle = {}

        #queues transformed by leviathon
        handler_bundle[Handler] = {"inqu":"main",
                                   "const":const, 
                                   "out_qu_dict":["merge","chaser"]}

        handler_bundle[parabam.merger.Handler] = {"inqu":"merge",
                                                  "const":const,
                                                  "out_qu_dict":[]}
        return handler_bundle

    def __get_output_paths__(self,input_paths,user_subsets,ensure_unique_output,**kwargs):
        output_paths = {}
        for salt_x,(input_path) in enumerate(input_paths)):
            output_paths[input_path] = {}
            for salt_y,subset in enumerate(user_subsets):
                output_paths[input_path][subset] = self.__get_path__(input_path,
                                                                     subset,
                                                                     ensure_unique_output,
                                                                     "%d%d" % (salt_x,salt_y))
        return output_paths

    def __get_path__(self,input_path,subset,ensure_unique_output,salt):
        head,ext = os.path.splitext(input_path)
        head,tail = os.path.split(head)
        unique = ""

        if ensure_unique_output:
            unique = "_%d%s" % (time.time(),salt)

        return os.path.join(".",self._temp_dir, "%s_%s%s%s" % (tail,subset,unique,ext))
            
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

    def __get_const_args__(self,**kwargs):
        args = super(Interface,self).__get_const_args__(**kwargs)
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