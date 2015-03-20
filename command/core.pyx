import time
import sys
import os
import gc
import shutil
import gzip

import pysam
import parabam

import Queue as Queue2

from itertools import izip
from multiprocessing import Queue,Process
from abc import ABCMeta, abstractmethod

#TODO: error handle incorrect return from user engine

class Task(parabam.core.Task):
    __metaclass__=ABCMeta

    def __init__(self,parent_bam,inqu,outqu,task_size,constants):
        super(Task, self).__init__(parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    task_size=task_size,
                                    constants=constants)

        self._engine = constants.user_engine
        self._user_constants = constants.user_constants

    @abstractmethod
    def __generate_results__(self,bam_iterator,**kwargs):
        #Must make call to self.__process_task_set__
        pass

    @abstractmethod
    def __post_run_routine__(self,**kwargs):
        pass

    @abstractmethod
    def __handle_engine_output__(self,engine_output,read):
        pass

    def __process_task_set__(self,bam_iterator):
        engine = self._engine
        next_read = bam_iterator.next 
        parent_bam = self._parent_bam
        cdef int size = self._task_size
        handle_output = self.__handle_engine_output__
        user_constants = self._user_constants
        counts = self._counts

        #StopIteration caught in parabam.core.Task.run
        for i in xrange(size):
            counts["total"] += 1
            read = next_read()
            engine_output = engine(read,user_constants,parent_bam)
            handle_output(engine_output,read)

    def __get_extension__(self,path):
        root,extension = os.path.splitext(path)
        return extension

class PairTask(Task):
    __metaclass__ = ABCMeta

    def __init__(self,parent_bam,inqu,outqu,task_size,constants):
        super(Task, self).__init__(parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    task_size=task_size,
                                    constants=constants)

        if constants.include_duplicates:
            self.__read_filter__ = self.__filter__
        else:
            self.__read_filter__ = self.__filter_duplicates__

    @abstractmethod
    def __generate_results__(self,**kwargs):
        #Must make call to self.__process_task_set__
        pass

    @abstractmethod
    def __handle_engine_output__(self,engine_output,read):
        pass

    def __process_task_set__(self,task_set):
        if type(task_set) == dict:
            self.__process_dict_task_set__(task_set)
        elif type(task_set) == list:
            self.__process_list_task_set__(task_set)

    def __process_dict_task_set__(self,task_set):
        user_constants = self._constants.user_constants
        engine = self.__engine__

        for qname,reads in task_set.items():
            output = engine( reads, user_constants, self._parent_bam)
            self.__handle_engine_output__(output, reads)

    def __process_list_task_set__(self,task_set):
        #speedup
        user_constants = self._constants.user_constants
        engine = self.__engine__
        query_loners = self.__query_loners__ 
        read_pair_decision = self.__handle_engine_output__
        #speedup

        loners = {}

        for read in task_set:
            read1,read2 = query_loners(read,loners)
            if read1:
                output = engine((read1,read2),user_constants,self._parent_bam)
                self.__handle_engine_output__(output,(read1,read2,))

        self.__stash_loners__(loners)
        del loners

    def __stash_loners__(self,loners):
        loner_path = self.__get_temp_path__("chaser")
        loner_file = pysam.AlignmentFile(loner_path,"wb",header=self._parent_bam.header)
        loner_count = 0
        for qname,read in loners.items():
            loner_count += 1
            loner_file.write(read)
        loner_file.close()
        self._system["chaser"] = (loner_count,loner_path,)

    def __query_loners__(self,read,loners):
        if self.__read_filter__(read):#Pair processing only handles primary pairs
            return None, None

        try:
            mate = loners[read.qname]
            del loners[read.qname]
            return read,mate

        except KeyError:
            loners[read.qname] = read
            return None,None

    def __filter__(self,read):
        return read.is_secondary or (read.flag & 2048 == 2048)

    def __filter_duplicates__(self,read):
        return read.is_secondary or (read.flag & 2048 == 2048) or read.is_duplicate

class Handler(parabam.core.Handler):
    __metaclass__=ABCMeta

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object constants,object pause_qu,dict out_qu_dict,object report=True):
        
        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                                     inqu=inqu,constants=constants,pause_qu=pause_qu,
                                     out_qu_dict=out_qu_dict)

        self._system_subsets = constants.system_subsets

        #Stage stores serve as a staging area
        #for packages to be sent on to subhandlers
        self._stage_stores = {} 
        for subset in self._system_subsets:
            self._stage_stores[subset] = []

        self._out_qu_dict = out_qu_dict
        self._stagecount = 0

    #Child classes MUST call super
    def __new_package_action__(self,new_package,**kwargs):
        results = new_package.results
        handle_dict = dict(results)
        #hack so as not record rescued paired reads twice in total
        if self._constants.pair_process and not new_package.parent_class == "Task":
            handle_dict["total"] = 0
        self.__auto_handle__(handle_dict,self._stats)

        if "system" in results.keys():
            self.__add_system_to_stage__(results["system"])
        
    def __add_system_to_stage__(self,system_results):
        for subset,(count,path) in system_results.items():
             self._stage_stores[subset].append((count,path))

    #Child classes MUST call super
    def __periodic_action__(self,iterations):
        self.__decide_if_finished__()

        for subset in self._system_subsets:
            if self.__test_stage_store__(subset):
                self.__add_staged_system_task__(results=self._stage_stores[subset],
                                        subset_type=subset)
                self.__update_stage_store__(subset)

    def __decide_if_finished__(self):
        if self._destroy:
            has_stage_tasks = False
            for subset,staged in self._stage_stores.items():
                if len(staged) > 0:
                    has_stage_tasks = True 
                    break
            if not has_stage_tasks:
                self._finished = True
                    
    def __update_stage_store__(self,subset):
        self._stagecount += 1
        self._stage_stores[subset] = []
        gc.collect()

    def __test_stage_store__(self,subset):
        if self._destroy:
            return True
        else:
            return len(self._stage_stores[subset]) > 20
    
    def __add_staged_system_task__(self,results,subset_type):
        if subset_type == "chaser":
            res = parabam.chaser.OriginPackage(results=results,
                                chaser_type="origin")
            self._out_qu_dict["chaser"].put(res)

    #Child classes MUST call super
    def __handler_exit__(self, **kwargs):
        if self._verbose:
            self.__standard_output__("\n[Status] Processing complete. %d reads processed" \
                                        % (self.__total_reads__(),))

class Interface(parabam.core.Interface):

    __metaclass__ = ABCMeta

    def __init__(self,temp_dir):
        super(Interface,self).__init__(temp_dir)

    def __get_module_and_vitals__(self,code_path):
        if os.getcwd not in sys.path:
            sys.path.append(os.getcwd())
        try:
            module = __import__(code_path, fromlist=[''])

        except ImportError:
            sys.stderr.write("[Error] parabam can't find user specified instructions"\
                             "\tEnsure instruction code is in current working directory\n")
            raise SystemExit

        user_engine = module.engine
        user_constants = {}
        if hasattr(module,"set_constants"):
            module.set_constants(user_constants)

        return module,user_engine,user_constants

    @abstractmethod
    def __get_handler_bundle__(self,queues,object constants,task_class,**kwargs):
        pass

    @abstractmethod
    def __get_destroy_handler_order__(self):
        '''Should return a list of hanlder classes in the 
        order that they are to be destroyed by the leviathon
        all hanlder class should have a corresponding bundle entry'''
        pass 

    @abstractmethod
    def __get_output_paths__(self,input_paths,**kwargs):
        pass

    @abstractmethod
    def __get_task_class__(self,**kwargs):
        pass

    @abstractmethod
    def __get_queue_names__(self,**kwargs):
        '''Predefine the names of queues in this implementation
        these names are converted to queue objects by the leviathon
        use these names when reference queus in processor and handler
        bundles'''
        pass

    @abstractmethod
    def __get_extra_const_args__(self,**kwargs):
        '''return empty dict if not needed'''
        pass

    def __get_const_args__(self,**kwargs):
        '''Any function overwriting this one MUST make a super call to this function.
        i.e args = super(<class_name>,self).__get_const_args__(**kwargs)'''

        args = {}
        args["temp_dir"] = self._temp_dir


        defaults = ["total_procs","task_size","verbose","user_constants",
                    "reader_n","user_constants", "user_engine", "fetch_region", 
                    "pair_process", "include_duplicates","input_is_sam","debug"]

        for key,val in kwargs.items():
            if key in defaults:
                args[key] = val

        if args["pair_process"]:
            args["system_subsets"] = ["chaser"]
        else:
            args["system_subsets"] = []

        return args

    def __update_final_output_paths__(self,input_path,output_paths,final_output_paths):
        if "global" in output_paths.keys():
            final_output_paths["global"].extend(output_paths["global"])

        final_output_paths[input_path] = []
        if input_path in output_paths.keys():
            output_path = output_paths[input_path]
            if type(output_path) == dict:
                for subset,path in output_path.items():
                    final_output_paths[input_path].append(path)    
            elif type(output_paths) == list:
                final_output_paths[input_path] = output_path
            elif type(output_paths) == str:
                final_output_paths[input_path].append(output_path)
            else:
                sys.stderr.write("[Error] Unexpected output path type: %s \n" \
                    % (type(output_paths[input_path])))
                raise SystemExit

    def __output_files_to_cwd__(self,final_output_paths):
        revised_output_paths = {}
        for master_path,out_paths in final_output_paths.items():
            revised_output_paths[ master_path ] = []
            for path in out_paths:
                head,tail = os.path.split(path)
                new_path = os.path.abspath(os.path.join(".",tail))
                real_path = self.__move_output_file__(path,new_path)
                revised_output_paths[master_path].append(real_path)
        return revised_output_paths

    def __move_output_file__(self,current_path,new_path):
        try:
            shutil.move(current_path,new_path)
            return new_path
        except shutil.Error,e: 
            root,ext = os.path.splitext( new_path )
            alternate_path = root + "_%d" % time.time() + ext

            sys.stderr.write("[Warning] Something went wrong when copying output to working directory.\n")
            sys.stdout.write("[Update]Trying to create output using unique filename:\n")
            sys.stdout.write("\t\t%s\n" % alternate_path)

            shutil.move(current_path,alternate_path)
            return alternate_path

    def __report_file_names__(self,final_output_paths,input_paths):
        sys.stdout.write("\n[Status] This run will output the following files:\n")
        for master_path,child_paths in final_output_paths.items():
            if master_path in input_paths:
                for path in child_paths:
                    root,name = os.path.split(path)
                    sys.stdout.write("\t%s\n" % (name,))            
        sys.stdout.write("\n")

    def run(self,**kwargs):

        input_paths = kwargs["input_bams"]

        const_args = self.__get_const_args__(**kwargs)
        const_args.update(self.__get_extra_const_args__(**kwargs))
        constants = parabam.core.Constants(**const_args)

        task_class = self.__get_task_class__(**kwargs)
        queue_names = self.__get_queue_names__(**kwargs)

        handler_order = self.__get_destroy_handler_order__()
        handler_bundle = self.__get_handler_bundle__(constants=constants,
                                                     task_class=task_class,
                                                     **kwargs)

        if constants.pair_process:
            self.__prepare_for_pair_processing__(handler_bundle,handler_order,
                                                 queue_names,constants)
        
        update_interval = self.__get_update_interval__(constants.verbose)

        leviathon = parabam.core.Leviathon(constants,handler_bundle,
                                           handler_order,queue_names,
                                           update_interval,task_class)

        final_output_paths = {"global":[]}

        for input_path in input_paths:

            output_paths = self.__get_output_paths__(input_path=input_path,**kwargs)
            self.__update_final_output_paths__(input_path,output_paths,
                                               final_output_paths)
            if constants.verbose: 
                self.__report_file_names__(final_output_paths,input_path)

            leviathon.run(input_path,output_paths)

        if not kwargs["keep_in_temp"]:
            final_output_paths = self.__output_files_to_cwd__(final_output_paths)
        
        self.__remove_empty_entries__(final_output_paths)
        return final_output_paths

    def __get_update_interval__(self,verbose):
        if verbose == 1: 
            return 60
        else:
            return 1

    def __remove_empty_entries__(self,final_output_paths):
        for master_path,child_paths in final_output_paths.items():
            if len(child_paths) == 0:
                del final_output_paths[master_path]

    def __prepare_for_pair_processing__(self,handler_bundle,handler_order,queue_names,constants):        
        handler_bundle[parabam.chaser.Handler] = {"inqu":"chaser",
                                                  "constants":constants,
                                                  "out_qu_dict":["main"]}

        for handler_class,handler_args in handler_bundle.items():
            if issubclass(handler_class,Handler):
                handler_bundle[handler_class]["out_qu_dict"].append("chaser")
        queue_names.append("chaser")
        handler_order.insert(0,parabam.chaser.Handler)

    def default_parser(self):
        parser = super(Interface,self).default_parser()

        parser.add_argument('--instruc','-i',metavar='INSTRUCTION',required=True
            ,help='The instruction file, written in python, that we wish'\
            'to carry out on the input BAM.')
        parser.add_argument('--input','-b',metavar='INPUT', nargs='+',required=True
            ,help='The file(s) we wish to operate on. Multipe entries should be separated by a single space')
        parser.add_argument('--debug',action="store_true",default=False,
            help="Only the first 5million reads will be processed")
        parser.add_argument('--pair',action="store_true",default=False
            ,help="A pair processor is used instead of a conventional processor")
        parser.add_argument('-r','--region',type=str,metavar="REGION",nargs='?',default=None
            ,help="The subset process will be run only on reads from this region\n"\
            "Regions should be colon seperated as specified by samtools (eg \'chr1:1000,5000\')")
        parser.add_argument('-d',action="store_false",default=True,
            help="parabam will not process reads marked duplicate.")
        parser.add_argument('-v', choices=[0,1,2],default=0,type=int,
            help="Indicate the amount of information output by the program:\n"\
            "\t0: No output [Default]\n"\
            "\t1: Total Reads Processsed\n"\
            "\t2: Detailed output")

        return parser
        
    @abstractmethod
    def get_parser(self):
        pass

#...happily ever after
