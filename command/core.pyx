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
    def __init__(self,object task_set, object outqu,object curproc,
        object destroy,object parent_bam,object const,str parent_class):

        super(Task, self).__init__(task_set=task_set,
            outqu=outqu,
            curproc=curproc*len(const.sources),
            destroy=destroy,
            parent_bam = parent_bam,
            const=const,
            parent_class=parent_class)

    @abstractmethod
    def __generate_results__(self):
        #Must make call to self.__process_task_set__
        pass

    @abstractmethod
    def __handle_engine_output__(self,engine_output,read):
        pass

    def __process_task_set__(self,task_set):
        engine = self.__engine__
        handle_engine = self.__handle_engine_output__
        user_constants = self.const.user_constants

        for read in task_set:
            decision = engine(read,user_constants,self._parent_bam)
            handle_engine(decision,read)

    def __get_temp_object__(self,path,ext):
        if ext == ".gz" or ext == ".gzip":
            return gzip.open(path,"wb")
        elif ext == ".txt":
            return open(path,"w")
        else:
            return pysam.AlignmentFile(path,"wb",header=self._parent_bam.header)

    def __get_extension__(self,path):
        root,extension = os.path.splitext(path)
        return extension

    def __get_temp_path__(self,typ,ext=".bam"):
        #self.pid ensures that the temp names are unique.
        return "%s/%s_%s_%s_parabam_temp%s" %\
            (self._temp_dir,self._source,typ,self.pid,ext)

    def __engine__(self,read,user_constants,parent):
        return self.const.user_engine(read,user_constants,parent)

class PairTask(Task):
    __metaclass__ = ABCMeta

    def __init__(self,object task_set,object outqu,object curproc,
                 object destroy,object parent_bam,object const,str parent_class):
        
        super(PairTask, self).__init__(task_set=task_set,
                                       outqu=outqu,
                                       curproc=curproc,
                                       destroy=destroy,
                                       parent_bam = parent_bam,
                                       const=const,
                                       parent_class=parent_class)

        if const.include_duplicates:
            self.__read_filter__ = self.__filter__
        else:
            self.__read_filter__ = self.__filter_duplicates__

    @abstractmethod
    def __generate_results__(self):
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
        user_constants = self.const.user_constants
        engine = self.__engine__

        for qname,reads in task_set.items():
            output = engine( reads, user_constants, self._parent_bam)
            self.__handle_engine_output__(output, reads)

    def __process_list_task_set__(self,task_set):
        #speedup
        user_constants = self.const.user_constants
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
        return read.is_secondary

    def __filter_duplicates__(self,read):
        return read.is_secondary or read.is_duplicate

class Handler(parabam.core.Handler):
    __metaclass__=ABCMeta

    def __init__(self,object inqu,object const,
                 object destroy_limit,object out_qu_dict):
        
        super(Handler,self).__init__(inqu,const,destroy_limit=destroy_limit)

        self._sources = const.sources
        self._system_subsets = const.system_subsets

        #Stage stores keep serve as a staging area
        #for packages to be sent on to subhandlers
        self._stage_stores = {}
        for source in self._sources:
            self._stage_stores[source] = {} 
            for subset in self._system_subsets:
                self._stage_stores[source][subset] = []

        self._out_qu_dict = out_qu_dict
        self._proc_flag_sent = False
        self._stagecount = 0

    #Child classes MUST call super
    def __new_package_action__(self,new_package,**kwargs):
        results = new_package.results
        handle_dict = dict(results)
        source = results["source"]
        #hack so as not record rescued paired reads twice in total
        if self.const.pair_process and not new_package.parent_class == "PairProcessor":
            handle_dict["total"] = 0
        self.__auto_handle__(handle_dict,self._stats,source)

        if "system" in results.keys():
            self.__add_system_to_stage__(results["system"],source)
        
    def __add_system_to_stage__(self,system_results,source):
        for subset,(count,path) in system_results.items():
             self._stage_stores[source][subset].append((count,path))

    def __auto_handle__(self,results,stats,source):
        if source not in stats.keys():
            #if we autoupdate before receiving a package
            stats[source] = {}
        super(Handler,self).__auto_handle__(results,stats[source])

    #Child classes MUST call super
    def __periodic_action__(self,iterations):
        for source in self._sources:
            for subset in self._system_subsets:
                if self.__test_stage_store__(source,subset):
                    self.__add_staged_system_task__(results=self._stage_stores[source][subset],
                                            subset_type=subset,source=source)
                    self.__update_stage_store__(source,subset)
                    
    def __update_stage_store__(self,source,subset):
        self._stagecount += 1
        self._stage_stores[source][subset] = []
        gc.collect()

    def __test_stage_store__(self,source,subset):
        #Test whether there are any temporary bams to be merged
        if subset == "chaser" and not self.__is_processing__() and not self._proc_flag_sent:
            self._proc_flag_sent = True
            return True
        else:
            return len(self._stage_stores[source][subset]) > 0
    
    def __is_processing__(self):
        return self._destroy_count < (self._destroy_limit - 1)

    def __add_staged_system_task__(self,results,subset_type,source,destroy=False):
        if subset_type == "chaser":
            res = parabam.chaser.OriginPackage(results=results,destroy=destroy,
                                source=source,chaser_type="origin",
                                processing= self.__is_processing__())
            self._out_qu_dict["chaser"].put(res)

    def __total_reads__(self):
        return sum(map(lambda s : self._stats[s]["total"],self._sources))

    def __send_destroy__(self,source,subset):
        self.__add_staged_system_task__(results=self._stage_stores[source][subset],
                            subset_type=subset,source=source,destroy=True)
        self.__update_stage_store__(source,subset)

    #Child classes MUST call super
    def __handler_exit__(self, **kwargs):
        if self._verbose:
            self.__standard_output__("\n[Status] Processing complete. %d reads processed" \
                                        % (self.__total_reads__(),))
        
        for source in self._sources:
            for subset in self._system_subsets:
                self.__send_destroy__(source,subset)

class Processor(parabam.core.Processor):

    def __init__(self,str source,object outqu,object const,object TaskClass,
                 object task_args,object debug=False):

        self._source = source
        super(Processor,self).__init__(outqu,const,TaskClass,task_args,debug)

    def __get_parent_bam__(self, master_file_path):
        return parabam.core.ParentAlignmentFile(master_file_path[self._source],self.const.input_is_sam)

    def __get_master_bam__(self,master_file_path):
        if self.const.input_is_sam:
            return pysam.Samfile(self.const.master_file_path[self._source],"r")
        else:
            return pysam.Samfile(self.const.master_file_path[self._source],"rb")

    def __add_to_collection__(self,master,alig,collection):
        collection.append(alig)

    def __get_next_alig__(self,master_bam):
        if not self.const.fetch_region:
            for alig in master_bam.fetch(until_eof=True):
                yield alig
        else:
            for alig in master_bam.fetch(region=self.const.fetch_region):
                yield alig

    def __pre_processor__(self,master_file_path):
        pass

class PairProcessor(Processor):

    def __init__(self,str source,object outqu,object const,object TaskClass,
                 object task_args,object inqu,object debug=False):

        super(PairProcessor,self).__init__(source,outqu,const,TaskClass,task_args,debug)
        self._inqu = inqu

    def __query_pause_qu__(self):
        try:
            pause = self._inqu.get(False)
            if pause:
                while True:
                    try:
                        pause = self._inqu.get(False)
                        if not pause:
                            break
                    except Queue2.Empty:
                        time.sleep(3)
        except Queue2.Empty:
            pass

        last = False   
        while True:
            try:
                pause = self._inqu.get(False)
                last = pause
            except Queue2.Empty:
                break
        return last

    def __start_task__(self,collection,destroy=False):
        while True:
            pause = self.__query_pause_qu__()
            if not pause:
                break
        super(PairProcessor,self).__start_task__(collection,destroy)

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
    def __get_processor_bundle__(self,queues,object const,task_class,**kwargs):
        pass

    @abstractmethod
    def __get_handler_bundle__(self,queues,object const,task_class,**kwargs):
        pass

    @abstractmethod
    def __get_master_file_path__(self,sources,input_paths,**kwargs):
        pass

    @abstractmethod
    def __get_output_paths__(self,sources,input_paths,**kwargs):
        pass

    @abstractmethod
    def __get_task_class__(self,**kwargs):
        pass

    @abstractmethod
    def __get_queues__(self,object const):
        '''predefine queue objects for communication 
        between processors and handlers. Should return
        a dict of queues'''
        pass

    def __get_const_args__(self,**kwargs):
        '''Any function overwriting this one MUST make a super call to this function.
        i.e args = super(<class_name>,self).__get_const_args__(**kwargs)'''

        args = {}
        args["temp_dir"] = self._temp_dir
        defaults = ["chunk","verbose","user_constants",
                    "user_constants", "user_engine", "fetch_region", 
                    "pair_process", "include_duplicates","input_is_sam"]
                    #TODO: Handle proc division __getmaxproc__
        for key,val in kwargs.items():
            if key in defaults:
                args[key] = val

        args["proc"] = self.__get_max_proc__(kwargs["proc"],
                                           kwargs["side_by_side"],
                                           kwargs["pair_process"])

        if args["pair_process"]:
            args["system_subsets"] = ["chaser"]
        else:
            args["system_subsets"] = []

        return args

    def __create_handlers__(self,handler_bundle):
        handler_objects = []
        for handler_class,handler_arguments in handler_bundle.items():
            handler_objects.append(handler_class(**handler_arguments))
        return handler_objects

    def __create_processors__(self,processor_bundle):
        processor_objects = []
        processor_class = processor_bundle["class"]
        del processor_bundle["class"]

        for source,args in processor_bundle.items():
            processor_objects.append(processor_class(**args))
        return processor_objects

    def __update_final_output_paths__(self,sources,master_file_path,output_paths,final_output_paths):
        if "global" in output_paths.keys():
            final_output_paths["global"].extend(output_paths["global"])

        for source in sources:
            master_path = master_file_path[source]
            final_output_paths[master_path] = []
            if source in output_paths.keys():
                output_path = output_paths[source]
                if type(output_path) == dict:
                    for subset,path in output_path.items():
                        final_output_paths[master_path].append(path)    
                elif type(output_paths) == list:
                    final_output_paths[master_path] = output_path
                elif type(output_paths) == str:
                    final_output_paths[master_path].append(output_path)
                else:
                    sys.stderr.write("[Error] Unexpected output path type: %s \n" \
                        % (type(output_paths[source])))
                    raise SystemExit
            else:
                continue

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

    def __add_run_details_to_const__(self,object const,master_file_path,output_paths,sources):
        #put the run sensistive information into const
        const.add("master_file_path",master_file_path)
        const.add("sources",sources)
        const.add("output_paths",output_paths)

    def run(self,**kwargs):

        input_bams = kwargs["input_bams"]
        side_by_side = kwargs["side_by_side"]

        const_args = self.__get_const_args__(**kwargs)
        const = parabam.core.Const(**const_args)
        final_output_paths = {"global":[]}
        task_class = self.__get_task_class__(**kwargs)

        for input_paths,sources in self.__get_group__(input_bams,multi=side_by_side):

            master_file_path = self.__get_master_file_path__(sources=sources,
                                                               input_paths=input_paths,**kwargs)
            output_paths = self.__get_output_paths__(sources=sources,
                                                     input_paths=input_paths,**kwargs)

            self.__update_final_output_paths__(sources,master_file_path,output_paths,final_output_paths)
            self.__add_run_details_to_const__(const=const,master_file_path=master_file_path,
                                              output_paths=output_paths,sources=sources)

            if const.verbose: 
                self.__report_file_names__(final_output_paths,input_paths)

            queues = self.__get_queues__(const,**kwargs)

            handler_bundle = self.__get_handler_bundle__(queues=queues,
                                                         const=const,
                                                         task_class=task_class,
                                                         **kwargs)

            processor_bundle = self.__get_processor_bundle__(queues=queues,
                                                             const=const,
                                                             task_class=task_class,
                                                             **kwargs)

            if const.pair_process:
                self.__prepare_for_pair_processing__(processor_bundle,
                                                     handler_bundle,
                                                     task_class,
                                                     const)

            handler_objects = self.__create_handlers__(handler_bundle)
            processor_objects = self.__create_processors__(processor_bundle)

            if const.verbose == 1: 
                update_interval = 199
            else:
                update_interval = 3
            lev = parabam.core.Leviathon(processor_objects,handler_objects,update_interval)
            lev.run()
            del lev

        if not kwargs["keep_in_temp"]:
            final_output_paths = self.__output_files_to_cwd__(final_output_paths)
        
        self.__remove_empty_entries__(final_output_paths)
        return final_output_paths

    def __remove_empty_entries__(self,final_output_paths):
        for master_path,child_paths in final_output_paths.items():
            if len(child_paths) == 0:
                del final_output_paths[master_path]

    def __prepare_for_pair_processing__(self,processor_bundle,handler_bundle,task_class,object const):

        chaser_qu = Queue()
        pause_qus = []
        main_qu = None

        for source in const.sources:
            newqu = Queue()
            processor_bundle[source]["inqu"] = newqu
            pause_qus.append(newqu)
        processor_bundle["class"] = PairProcessor

        for entry_class,args in handler_bundle.items():
            if issubclass(entry_class,Handler):
                args["destroy_limit"] += 1
                main_qu = args["inqu"]
                args["out_qu_dict"]["chaser"] = chaser_qu
        
        handler_bundle[parabam.chaser.Handler] = {"inqu":chaser_qu,
                                                  "pause_qus":pause_qus,
                                                  "mainqu":main_qu,
                                                  "const":const,
                                                  "destroy_limit":1,
                                                  "TaskClass":task_class}

    def __get_max_proc__(self,proc,input_size,pair_process):
        max_proc = (proc // input_size)
        if pair_process:
            max_proc = max_proc // 2
        if max_proc == 0:
            sys.stderr.write('[Error] Too many side-by-side inputs and too few processors.\n')
            sys.stderr.write('\tNot enough processors to go round. Reduce -s or increase -p')
            sys.stderr.flush()
            raise SystemExit()
        return max_proc

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
        parser.add_argument('-s',type=int,metavar="INT",nargs='?',default=2
            ,help="Further parralise by running this many samples simultaneously [Default 2]")
        parser.add_argument('-d',action="store_true",default=False,
            help="parabam will not process reads marked PCR duplicate.")
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