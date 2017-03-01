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

#TODO: error handle incorrect return from user rule

class Task(parabam.core.Task):
    __metaclass__=ABCMeta

    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants):
        super(Task, self).__init__(parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    statusqu=statusqu,
                                    task_size=task_size,
                                    constants=constants)

        self._user_rule = constants.rule
        self._user_constants = constants.constants

    @abstractmethod
    def __handle_rule_output__(self,rule_output,read):
        pass

    def __pre_run_routine__(self,iterator,**kwargs):
        self._system = {}

    def __post_run_routine__(self,**kwargs):
        pass

    def __process_task_set__(self,iterator):
        user_rule = self._user_rule
        next_read = iterator.next 
        parent_bam = self._parent_bam
        handle_output = self.__handle_rule_output__
        user_constants = self._user_constants
        
        #StopIteration caught in parabam.core.Task.run
        for i in xrange(self._task_size):    
            read = next_read()
            rule_output = user_rule(read,user_constants,parent_bam)
            handle_output(rule_output,read)

    def __get_extension__(self,path):
        root,extension = os.path.splitext(path)
        return extension

class PairTask(Task):
    __metaclass__ = ABCMeta

    def __init__(self,
                 parent_bam,
                 inqu,
                 outqu,
                 statusqu,
                 task_size,
                 constants):
        super(PairTask, self).__init__(parent_bam=parent_bam,
                                       inqu=inqu,
                                       outqu=outqu,
                                       statusqu=statusqu,
                                       task_size=task_size,
                                       constants=constants)

        if constants.include_duplicates:
            self.__read_filter__ = self.__filter__
        else:
            self.__read_filter__ = self.__filter_duplicates__

    @abstractmethod
    def __handle_rule_output__(self,rule_output,read):
        pass

    def __pre_run_routine__(self,iterator):
        super(PairTask,self).__pre_run_routine__(iterator)
        self._loners = {}

    def __post_run_routine__(self,**kwargs):
        super(PairTask,self).__post_run_routine__()
        self.__stash_loners__(self._loners)
        del self._loners

    def __process_task_set__(self,iterator):
        
        next_read = iterator.next 
        query_loners = self.__query_loners__ 
        handle_output = self.__handle_rule_output__

        user_rule = self._user_rule
        user_constants = self._user_constants
        parent_bam = self._parent_bam

        counts = self._counts
        loners = self._loners

        cdef int size = self._task_size

        #StopIteration caught in parabam.core.Task.run
        for i in xrange(size):
            read = next_read()
            read1,read2 = query_loners(read,loners)

            if read1:
                rule_output = user_rule((read1,read2), 
                                         user_constants,
                                         parent_bam)
                handle_output(rule_output, (read1, read2,))

    def __stash_loners__(self,loners):
        loner_count = 0
        loner_path = self.__get_temp_path__("chaser")
        loner_file = pysam.AlignmentFile(loner_path,"wb",
                                         header=self._parent_bam.header)

        for qname,read in loners.items():
            loner_count += 1
            loner_file.write(read)
        self._system["chaser"] = (loner_count,loner_path,)
        loner_file.close()

    def __query_loners__(self,read,loners):
        #Pair processing only handles primary pairs
        if self.__read_filter__(read):
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
        return read.is_secondary or \
                (read.flag & 2048 == 2048) \
                    or read.is_duplicate

class ByCoordTask(Task):
    __metaclass__ = ABCMeta

    def __init__(self,
                 parent_bam,
                 inqu,
                 outqu,
                 statusqu,
                 task_size,
                 constants):
        super(ByCoordTask, self).__init__(parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    statusqu=statusqu,
                                    task_size=task_size,
                                    constants=constants)

        #TODO: Secondary alignments???
        self.__read_filter__ = self.__filter__

    @abstractmethod
    def __handle_rule_output__(self,rule_output,read):
        pass

    def __process_task_set__(self,iterator):
        user_rule = self._rule
        next_read = iterator.next 
        parent_bam = self._parent_bam

        handle_output = self.__handle_rule_output__
        user_constants = self._user_constants

        position_max = 10000

        task_pos = None
        reads = []

        self._task_size = 0
        current_positions = 0

        #StopIteration caught in parabam.core.Task.run
        while current_positions <= position_max:
            read = next_read()
            if task_pos is None:
                task_pos = read.pos

            if task_pos == read.pos:
                reads.append(read)
                
            else:
                rule_output = user_rule(reads,user_constants,parent_bam)        
                handle_output(rule_output,reads)
                self._task_size += len(reads)
                
                del reads
                reads = [read]
                task_pos = read.pos
                current_positions += 1

    def __filter__(self,read):
        return read.is_secondary or (read.flag & 2048 == 2048)

class Handler(parabam.core.Handler):
    __metaclass__=ABCMeta

    def __init__(self, 
                object parent_bam, 
                object output_paths,
                object inqu,
                object constants,
                object pause_qus,
                dict out_qu_dict,
                object report=True):
        
        super(Handler,self).__init__(parent_bam=parent_bam,
                                     output_paths=output_paths,
                                     inqu=inqu,
                                     constants=constants,
                                     pause_qus=pause_qus,
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
    def __new_package_action__(self, new_package, **kwargs):
        results = new_package.results
        self.__auto_handle__(results,self._stats)

        if "system" in results.keys():
            self.__add_system_to_stage__(results["system"])
        
    def __add_system_to_stage__(self, system_results):
        for subset,(count,path) in system_results.items():
             self._stage_stores[subset].append((count,path))

    #Child classes MUST call super
    def __periodic_action__(self, iterations):
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
                    
    def __update_stage_store__(self, subset):
        self._stagecount += 1
        self._stage_stores[subset] = []
        gc.collect()

    def __test_stage_store__(self, subset):
        if (not self._processing) or self._destroy:
            return len(self._stage_stores[subset]) > 0
        else:
            return len(self._stage_stores[subset]) > 10
    
    def __add_staged_system_task__(self, results, subset_type):
        if subset_type == "chaser":
            res = parabam.chaser.ChaserResults(results=results,
                                chaser_type="origin")
            self._out_qu_dict["chaser"].put(res)

    #Child classes MUST call super
    def __handler_exit__(self, **kwargs):
        if self._verbose:
            self.__standard_output__("\t- Total reads processed: %d" \
                                        % (self.__total_reads__(),))

class ByCoordFileReader(parabam.core.FileReader):

    def __init__(self,str input_path,int proc_id,object outqu,
                 int task_n,object constants,object Task,
                 object pause_qu,object inqu = None):

        super(ByCoordFileReader,self).__init__(input_path=input_path,
                                                proc_id=proc_id,
                                                outqu=outqu,
                                                task_n=task_n,
                                                constants=constants,
                                                Task=Task,
                                                pause_qu=pause_qu,
                                                inqu=inqu,)
    def __bam_generator__(self,bam_file):
        cdef int proc_id = self._proc_id
        cdef int reader_n = self._reader_n
        cdef int iterations = 0
        cdef int task_size = self._task_size

        cdef int position_max = 10000
        
        parent_iter = bam_file.fetch(until_eof = True)

        one_ahead_bam = pysam.AlignmentFile(bam_file.filename,"rb")
        one_ahead_iter = one_ahead_bam.fetch(until_eof = True)
        one_ahead_iter.next()

        while True:
            try:
                if iterations % reader_n == proc_id:
                    yield iterations
                
                observed_positions = 0
                while True:

                    current_read_pos = parent_iter.next().pos
                    try:
                        one_ahead_pos = one_ahead_iter.next().pos                        
                    except StopIteration:
                        one_ahead_pos = -1

                    if not current_read_pos == one_ahead_pos:
                        observed_positions += 1
                        if observed_positions == position_max:
                            break

                iterations += 1

            except StopIteration:
                break
        return

    def __get_parent_iter__(self,parent_bam):
        return parent_bam

class Interface(parabam.core.Interface):

    __metaclass__ = ABCMeta

    def __init__(self,pair_process = False,
                      coord_process = False,
                      include_duplicates = True,
                      debug = False,
                      ensure_unique_output=False,
                      **kwargs):

        super(Interface,self).__init__(**kwargs)

        if self.cmd_run == False:
            self.pair_process = pair_process

            self.coord_process = False
            if not self.pair_process:
                self.coord_process = coord_process

            self.include_duplicates = include_duplicates
            self.debug = debug
            self.ensure_unique_output = ensure_unique_output

    def __get_module_and_vitals__(self,code_path):
        if os.getcwd not in sys.path:
            sys.path.append(os.getcwd())

        code_path = code_path.replace(".py","")

        # TODO: I've seen an error where a bug in the imported package
        #       causes parabam to throw this message. Needs further
        #       exploration
        try:
            module = __import__(code_path, fromlist=[''])
        except ImportError: 
            sys.stderr.write("[Error] parabam can't find user specified instructions\n"\
                              "\tEnsure instruction code is in current working directory\n")
            raise SystemExit

        rule = module.rule
        constants = {}
        if hasattr(module,"set_constants"):
            module.set_constants(constants)

        return module,rule,constants

    def __cmd_args_to_class_vars__(self):

        super(Interface,self).__cmd_args_to_class_vars__()
        
        self.pair_process = self.cmd_args.pair
        self.coord_process = False

        if not self.pair_process:
            self.coord_process = self.cmd_args.coord

        self.include_duplicates = self.cmd_args.d
        self.debug = self.cmd_args.debug

    @abstractmethod
    def __get_handler_bundle__(self,**kwargs):
        pass

    @abstractmethod
    def __get_destroy_handler_order__(self):
        '''Should return a list of hanlder classes in the 
        order that they are to be destroyed by the leviathon
        all hanlder class should have a corresponding bundle entry'''
        pass 

    @abstractmethod
    def __get_output_paths__(self,input_path,final_output_paths,**kwargs):
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
    def __instalise_final_output__(self,**kwargs):
        '''Allow sub-classes to instalise the dictionary
        of output paths that are returned to the caller
        of the run function.
        '''
        pass

    def __get_const_args__(self,**kwargs):
        '''Any function overwriting this one MUST make a 
            super call to this function.

            i.e args = super(<class_name>,self).__get_const_args__(**kwargs)'''

        args = {}
        
        #Dump all non function variables from this interface into constants
        for attribute in dir(self):
            is_function = hasattr(getattr(self,attribute),"__call__")
            if not is_function and not attribute.startswith("_"):
                args[attribute] = getattr(self,attribute)

        for key,val in kwargs.items():
            args[key] = val

        if self.pair_process:
            args["system_subsets"] = ["chaser"]
        else:
            args["system_subsets"] = []

        return args

    def __update_final_output_paths__(self,
                                       input_path,
                                       output_paths,
                                       final_output_paths):
        final_output_paths.update(output_paths)
    
    def __output_files_to_cwd__(self,final_output_paths):
        revised_output_paths = {}
        for input_path, analyses in final_output_paths.items():
            revised_output_paths[input_path] = {}
            for analysis_name, analysis_path in analyses.items():
                head,tail = os.path.split(analysis_path)
                new_path = os.path.abspath(os.path.join(".",tail))

                real_path = self.__move_output_file__(analysis_path,new_path)
                revised_output_paths[input_path][analysis_name] = real_path

        return revised_output_paths

    def __move_output_file__(self,current_path,new_path):
        try:
            shutil.move(current_path,new_path)
            return new_path
        except shutil.Error,e: 
            root,ext = os.path.splitext( new_path )
            alternate_path = root + "_%d" % time.time() + ext

            sys.stderr.write("WARNING: Something went wrong when copying output to working directory.\n")
            sys.stdout.write("... Trying to create output using unique filename:\n")
            sys.stdout.write("\t\t%s\n" % alternate_path)

            shutil.move(current_path,alternate_path)
            return alternate_path

    def __report_file_names__(self,final_output_paths,input_paths):
        sys.stdout.write("\t- This run will output the following files:\n")
        for master_path,analyses in final_output_paths.items():
            for name,path in analyses.items():
                root,name = os.path.split(path)
                sys.stdout.write("\t\t+ %s\n" % (name,))            

    def run(self,input_paths,**kwargs):
        self.__temp_dir_instalise__()
        self.__introduce__()

        const_args = self.__get_const_args__(**kwargs)
        constants = parabam.core.Constants(**const_args)

        task_class = self.__get_task_class__(**kwargs)
        queue_names = self.__get_queue_names__(**kwargs)

        handler_order = self.__get_destroy_handler_order__()
        handler_bundle = self.__get_handler_bundle__(**kwargs)

        if constants.pair_process:
            self.__prepare_for_pair_processing__(handler_bundle,
                                                 handler_order,
                                                 queue_names,
                                                 constants,
                                                 task_class)
        
        update_interval = self.__get_update_interval__(constants.verbose)

        filereader_class = self.__get_filereader_class__(**kwargs)

        leviathon = parabam.core.Leviathan(constants,handler_bundle,
                                           handler_order,queue_names,
                                           update_interval,task_class,
                                           filereader_class)

        final_output_paths = self.__instalise_final_output__(**kwargs)

        for input_path in input_paths:
            output_paths = self.__get_output_paths__(input_path,
                                                     final_output_paths,
                                                     **kwargs)

            self.__update_final_output_paths__(input_path,output_paths,
                                               final_output_paths)

            if constants.verbose and not self.keep_in_temp: 
                self.__report_file_names__(final_output_paths,input_path)

            leviathon.run(input_path,output_paths)

        if not self.keep_in_temp:
            final_output_paths = self.__output_files_to_cwd__(final_output_paths)
        
        self.__remove_empty_entries__(final_output_paths)

        self.__goodbye__()
        self.interface_exit()

        return final_output_paths

    def __get_update_interval__(self,verbose):
        if verbose == 1: 
            return 8000
        else:
            return 200

    def __get_filereader_class__(self,**kwargs):
        if self.coord_process:
            return ByCoordFileReader
        else:
            return parabam.core.FileReader

    def __remove_empty_entries__(self,final_output_paths):
        for master_path,child_paths in final_output_paths.items():
            if len(child_paths) == 0:
                del final_output_paths[master_path]

    def __prepare_for_pair_processing__(self,handler_bundle,handler_order,
                                            queue_names,constants,Task):        
        handler_bundle[parabam.chaser.Handler] = {"inqu":"chaser",
                                                  "constants":constants,
                                                  "out_qu_dict":["main"],
                                                  "Task":Task}

        for handler_class,handler_args in handler_bundle.items():
            if issubclass(handler_class,Handler):
                handler_bundle[handler_class]["out_qu_dict"].append("chaser")
        queue_names.append("chaser")
        handler_order.insert(0,parabam.chaser.Handler)
        constants.total_procs = constants.total_procs / 2

    def default_parser(self):
        parser = super(Interface,self).default_parser()

        parser.add_argument('--rule','-i',metavar='RULE',required=True
            ,help='The file containing the rule, written in python,\n'\
            'that we wish to apply to the input BAM.')
        parser.add_argument('--input','-b',metavar='INPUT', 
                            nargs='+',required=True
            ,help='The file(s) we wish to operate on.\n'\
            ' Multiple entries should be separated by a single space')
        parser.add_argument('--debug',action="store_true",default=False,
            help="Only the first 5million reads will be processed")
        parser.add_argument('--pair',action="store_true",default=False
            ,help="A pair processor is used instead of a conventional processor")
        parser.add_argument('--coord',action="store_true",default=False
            ,help="Engines recieve a list of reads which all map to"\
                  " the same starting position.\n This mode is not"\
                  " compatible with the `--pair` option")
        parser.add_argument('-r','--region',type=str,
                                            metavar="REGION",
                                            nargs='?',
                                            default=None
            ,help="The process will be run only on reads from\n"\
            "this region. Regions should be colon separated as\n"\
            "specified by samtools (eg \'chr1:1000,5000\')")
        parser.add_argument('-d',action="store_false",default=True,
            help="parabam will not process reads marked duplicate.")

        return parser
        
    @abstractmethod
    def get_parser(self):
        pass

#...happily ever after
