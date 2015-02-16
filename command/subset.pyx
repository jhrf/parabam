import pysam
import pdb
import parabam
import time
import sys
import os
import gc
import shutil
import gzip

import Queue as Queue2

from itertools import izip
from collections import deque, Counter
from multiprocessing import Queue,Process
from abc import ABCMeta, abstractmethod
from parabam.support import HandlerMerge,MergePackage
from parabam.chaser import HandlerChaser,OriginPackage

class TaskSubset(parabam.core.Task):

    def __init__(self, object task_set, object outqu, object curproc,
                 object destroy, object const, str parent_class, str source):
        super(TaskSubset, self).__init__(task_set=task_set,
                                         outqu=outqu,
                                         curproc=curproc*len(const.sources),
                                         destroy=destroy,
                                         const=const,
                                         parent_class=parent_class)
        self._source = source
        self._master_file_path = const.master_file_path[self._source]
        self._subset_types = const.subset_types
        self._counts = {}
        self._temp_paths = {}
        self._temp_objects = {}

    def __generate_results__(self):
        master = pysam.Samfile(self._master_file_path,"rb")

        cdef dict temp_paths = self._temp_paths
        cdef dict temp_objects = self._temp_objects
        cdef dict counts = self._counts

        for subset in self._subset_types:
            ext = self.__get_extension__(self.const.output_paths[self._source][subset])
            temp_paths[subset] = self.__get_temp_path__(subset,ext)
            temp_objects[subset] = self.__get_temp_object__(temp_paths[subset],ext,master)
            counts[subset] = 0

        self.__handle_task_set__(self._task_set,master)

        #Close the master bamfile
        for subset,file_object in temp_objects.items():
            file_object.close() #close all the other bams
        del self._temp_objects
        master.close()

        results = {}
        results["source"] = self._source
        results["temp_paths"] = temp_paths
        results["counts"] = counts

        return results

    def __get_temp_object__(self,path,ext,master):
        if ext == ".gz" or ext == ".gzip":
            return gzip.open(path,"wb")
        elif ext == ".txt":
            return open(path,"w")
        else:
            return pysam.Samfile(path,"wb",template=master)

    def __get_extension__(self,path):
        root,extension = os.path.splitext(path)
        return extension

    def __write_to_subset_bam__(self,subset_type,read):
        self._counts[subset_type] += 1
        self._temp_objects[subset_type].write(read)

    def __handle_task_set__(self,task_set,master):

        engine = self.__engine__
        user_constants = self.const.user_constants
        subset_types = self.const.subset_types

        subset_write = self.__write_to_subset_bam__

        for read in task_set:
            subset_decision = engine(read,user_constants,master)
            
            if type(subset_decision) == bool:
                if subset_decision:
                    subset_write(subset_types[0],read)          
            elif type(subset_decision) == list:
                for subset,cur_read in subset_decision:
                    subset_write(subset,cur_read)
            elif type(subset_decision) == int:
                if not subset_decision == -1:
                    subset_write(subset_types[subset_decision],read)
            elif type(subset_decision) == None:
                pass
            else:
                sys.stdout.write("[ERROR] Unrecognised return type from user engine!\n")
                print subset_decision
                print type(subset_decision)
                print self.__class__
                sys.stdout.flush()

    def __get_temp_path__(self,typ,ext=".bam"):
        #self.pid ensures that the temp names are unique.
        return "%s/%s_%s_%s_parabam_temp%s" % (self._temp_dir,self._source,typ,self.pid,ext)

    def __engine__(self,read,user_constants,master):
        return self.const.user_engine(read,user_constants,master)

class PairTaskSubset(TaskSubset):
    def __init__(self,object task_set,object outqu,object curproc,
                 object destroy,object const,str parent_class,str source):
        
        super(PairTaskSubset, self).__init__(task_set=task_set,
                                             outqu=outqu,
                                             curproc=curproc,
                                             destroy=destroy,
                                             const=const,
                                             parent_class=parent_class,
                                             source=source)
        if const.include_duplicates:
            self.__read_filter__ = self.__filter__
        else:
            self.__read_filter__ = self.__filter_duplicates__

    def __handle_task_set__(self,task_set,master):
        engine = self.__engine__
        user_constants = self.const.user_constants

        loners = {}

        if type(task_set) == dict:
            self.__handle_dict_task_set__(task_set,engine,user_constants,master)

        elif type(task_set) == list:
            self.__handle_list_task_set__(task_set,engine,user_constants,master)

    def __handle_dict_task_set__(self,task_set,engine,user_constants,master):
        for qname,(read1,read2) in task_set.items():
            self.__read_pair_decision__(read1,read2,engine,user_constants,master)

    def __handle_list_task_set__(self,task_set,engine,user_constants,master):
        user_constants = self.const.user_constants
        subset_types = self.const.subset_types

        #speedup
        query_loners = self.__query_loners__ 
        read_pair_decision = self.__read_pair_decision__
        #speedup

        loners = {}

        for read in task_set:
            read1,read2 = query_loners(read,loners)
            if read1:
                read_pair_decision(read1,read2,engine,user_constants,master)

        self.__stash_loners__(loners)
        del loners

    def __stash_loners__(self,loners):
        subset_write = self.__write_to_subset_bam__
        for qname,read in loners.items():
            subset_write("index",read)

    def __read_pair_decision__(self,read1,read2 ,engine,user_constants,master):
        subset_decision = engine((read1,read2),user_constants,master)
        if type(subset_decision) == list:
            for subset,cur_read in subset_decision:
                self.__write_to_subset_bam__(subset,cur_read)
        else:
            sys.stdout.write("[ERROR] Unrecognised return type from user engine!\n")
            sys.stdout.write("[ERROR] In paired subset mode engine must return like-so:\n")
            sys.stdout.write("\t[ (subset_type,read), ... ]")
            sys.stdout.flush()

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

class HandlerSubset(parabam.core.Handler):

    def __init__(self,object inqu,object outqu,object const,
                 object destroy_limit,object chaserqu):
        
        super(HandlerSubset,self).__init__(inqu,const,destroy_limit=destroy_limit)

        self._sources = const.sources
        self._subset_types = const.subset_types

        #Setup stores and connection to merge proc
        self._merge_stores = {}
        for source in self._sources:
            self._merge_stores[source] = {} 
            for subset in self._subset_types:
                self._merge_stores[source][subset] = []

        self._mergequeue = outqu
        self._mergecount = 0

        self._proc_flag_sent = False
        self._chasequeue = chaserqu

    def __new_package_action__(self,new_package,**kwargs):
        results = new_package.results
        handle_dict = dict(results)
        source = results["source"]

        #hack so as not record rescued paired reads twice in total
        if self.const.pair_process and not new_package.parent_class == "ProcessorSubsetPair":
            handle_dict["total"] = 0
        if "index" in results["counts"].keys():
            handle_dict["counts"] = dict(handle_dict["counts"])
            del handle_dict["counts"]["index"]

        self.__auto_handle__(handle_dict,self._stats,source)       
        for subset in self._subset_types:
            self._merge_stores[source][subset].append((results["counts"][subset],results["temp_paths"][subset],))

    def __auto_handle__(self,results,stats,source):
        if source not in stats.keys():
            stats[source] = {}
        super(HandlerSubset,self).__auto_handle__(results,stats[source])

    def __periodic_action__(self,iterations):
        for source in self._sources:
            for subset in self._subset_types:
                if self.__test_merge_store__(source,subset):

                    self.__add_merge_task__(name=self._output_paths[source][subset],
                                            results=self._merge_stores[source][subset],
                                            subset_type=subset,source=source)
                    self.__update_merge_store__(source,subset)
                    
    def __update_merge_store__(self,source,subset):
        self._mergecount += 1
        self._merge_stores[source][subset] = []
        gc.collect()

    def __test_merge_store__(self,source,subset):
        #Test whether there are any temporary bams to be merged
        if subset == "index" and not self.__is_processing__() and not self._proc_flag_sent:
            self._proc_flag_sent = True
            return True
        else:
            return len(self._merge_stores[source][subset]) > 0
    
    def __is_processing__(self):
        return self._destroy_count < (self._destroy_limit - 1)

    def __add_merge_task__(self,name,results,subset_type,source,destroy=False):
        if subset_type == "index":
            res = OriginPackage(results=results,destroy=destroy,
                                source=source,chaser_type="origin",
                                processing= self.__is_processing__()  ) 
            self._chasequeue.put(res)
        else:
            res = MergePackage(name=name,results=results,
                        subset_type=subset_type,source=source,
                        destroy=destroy)
            self._mergequeue.put(res)

    def __total_reads__(self):
        return sum(map(lambda s : self._stats[s]["total"],self._sources))

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

    def __handler_exit__(self, **kwargs):
        if self._verbose:
            self.__standard_output__("\n[Status] Processing complete. %d reads processed" % (self.__total_reads__(),))
            self.__standard_output__("[Status] Waiting for merge operation to finish")

        #Kill the handlers handler
        for source in self._sources:
            for subset in self._subset_types:
                self.__add_merge_task__(name=self._output_paths[source][subset],
                                results=self._merge_stores[source][subset],subset_type=subset,
                                source=source,destroy=True)
                self.__update_merge_store__(source,subset)

        self.__write_counts_csv__()

class ProcessorSubset(parabam.core.Processor):

    def __init__(self,object outqu,object const,object TaskClass,object task_args,object debug=False):
        # if const.fetch_region:
        #   debug = False 
        super(ProcessorSubset,self).__init__(outqu,const,TaskClass,task_args,debug)
        self._source = task_args[0] #Defined in the run function within Interface

    def __get_master_bam__(self,master_file_path):
        return pysam.Samfile(master_file_path[self._source],"rb")

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

class ProcessorSubsetPair(ProcessorSubset):

    def __init__(self,object outqu,object const,object TaskClass,object task_args,object inqu,object debug=False):
        super(ProcessorSubsetPair,self).__init__(outqu,const,TaskClass,task_args,debug)
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
        super(ProcessorSubsetPair,self).__start_task__(collection,destroy)

class Interface(parabam.core.UserInterface):
    """The interface to parabam subset.
    Users will primarily make use of the ``run`` function."""

    def __init__(self,temp_dir):
        super(Interface,self).__init__(temp_dir)

    def run_cmd(self,parser):
        cmd_args = parser.parse_args()

        verbose = cmd_args.v

        module,user_engine,user_constants = self.__get_module_and_vitals__(cmd_args.instruc)

        if hasattr(module,"get_subset_types"):
            subset_types = module.get_subset_types()
        else:
            subset_types = ["subset"]

        self.run(
            input_bams=cmd_args.input,
            proc= cmd_args.p,
            chunk= cmd_args.c,
            verbose= verbose,
            subset_types= subset_types,
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
    
    def run(self,input_bams,proc,chunk,subset_types,
            user_constants,user_engine,fetch_region=None,side_by_side=2,
            keep_in_temp=False,engine_is_class=False,verbose=False,
            pair_process=False,include_duplicates=False,debug=False,
            ensure_unique_output=False,output_counts=False):


        #AT SOME POINT WE SHOULD HANDLE UNSORTED BAMS. EITHER HERE OR AT THE PROCESSOR
        final_files = []

        outputs = self.__get_outputs__(input_bams,unique=ensure_unique_output) 

        if pair_process and not "index" in subset_types:
            subset_types.append("index")

        for input_group,output_group in self.__get_group__(input_bams,outputs,multi=side_by_side):
            
            output_paths = dict([(source,{}) for source in output_group])
            master_file_path = {}

            for master_path,source in zip(input_group,output_group):
                master_file_path[source] = master_path
                for subset_type in subset_types:
                    output_paths[source][subset_type] = "%s/%s_%s.bam" % (self._temp_dir,source,subset_type,)
                    
            if verbose: self.__report_file_names__(output_paths)

            proc_per_processor = self.__get_max_proc__(proc,len(input_group),pair_process)

            const = parabam.core.Const(output_paths=output_paths,
                                temp_dir=self._temp_dir,
                                master_file_path=master_file_path,
                                chunk=chunk,proc=proc_per_processor,
                                verbose=verbose,
                                subset_types=subset_types,
                                sources=output_group,
                                user_constants=user_constants,
                                user_engine=user_engine,
                                fetch_region=fetch_region,
                                pair_process=pair_process,
                                include_duplicates=include_duplicates,
                                output_counts=output_counts)

            task_qu = Queue()
            processors,task_class,pause_qus = self.__create_processors__(task_qu,const,debug,engine_is_class)
            handlers = self.__create_handlers__(task_qu,pause_qus,const,task_class,proc)

            if verbose == 1: 
                update_interval = 199
            else:
                update_interval = 3

            lev = parabam.core.Leviathon(processors,handlers,update_interval)
            lev.run()
            del lev

            #Move the complete BAMs etc out of the temp_dir to the working dir
            #Only do this if we custom generated the file locations.
            if keep_in_temp:
                for source,subset_paths in output_paths.items():
                    for subset,path in subset_paths.items():
                        if not subset == "index": 
                            final_files.append(path)
            else:
                final_files.extend(self.__move_output_files__(output_paths))
            
            gc.collect()

        return final_files

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

    def __get_outputs__(self,input_bam_paths,unique):
        outputs = []
        if unique:
            unique_str = "_%d%s" % (os.getpid(),str(time.time()).replace(".",""),)
        else:
            unique_str = ""
            
        return [ "%s%s" % (self.__get_basename__(b),unique_str) for b in input_bam_paths ]

    def __seperate_subset_and_index__(self,output_paths):
        index_paths = {}
        subset_paths = {}
        for source,subset_paths in output_paths.items():
            index_paths[source] = {}
            subset_paths[source] = {}
            for subset,path in subset_paths.items():
                if subset == "index":
                    index_paths[source][subset] = path
                else:
                    subset_paths[source][subset] = path
        return index_paths,subset_paths

    def __create_handlers__(self,task_qu,pause_qus,object const,task_class,proc):
        handlers = []
        merge_qu = Queue()
        chaser_qu = Queue()

        if const.pair_process:
            #This modifies the expected destroy_count
            #Main handler waits for chaser
            #Merge handler doesn't expect index destroy
            destroy_modifier = 1
        else:
            destroy_modifier = 0

        handlers.append(HandlerSubset(inqu=task_qu,outqu=merge_qu,
                        const=const,chaserqu=chaser_qu,
                        destroy_limit=len(const.sources) + destroy_modifier))

        handlers.append(HandlerMerge(inqu=merge_qu,const=const,
                        destroy_limit=len(const.sources)*(len(const.subset_types) - destroy_modifier) ))

        if const.pair_process:
            handlers.append(HandlerChaser(inqu=chaser_qu,pause_qus=pause_qus,mainqu=task_qu,
                                const=const,destroy_limit=1,TaskClass=task_class,chaser_task_max= (proc//2)))

        return handlers

    def __create_processors__(self,task_qu,object const,debug,engine_is_class):
        processors = []
        pause_qus = []

        if const.pair_process:
            task_class = PairTaskSubset
        else:
            task_class = TaskSubset

        if engine_is_class:
            if not const.pair_process and not issubclass(const.user_engine,TaskSubset):
                raise_exception = True
            elif const.pair_process and not issubclass(const.user_engine,PairTaskSubset):
                raise_exception = True
            else:
                raise_exception = False
                
            if raise_exception:
                raise Exception("[ERROR]\tUser engine class must inherit %s\n" \
                    % (task_class.__class__,))
            else:
                task_class = const.user_engine

        for source in const.sources:
            if const.pair_process:
                pause_qus.append(Queue())
                processors.append(ProcessorSubsetPair(outqu=task_qu,
                    const=const,
                    TaskClass=task_class,
                    task_args=[source],
                    inqu = pause_qus[-1],
                    debug = debug))
            else:
                processors.append(ProcessorSubset(outqu=task_qu,
                    const=const,
                    TaskClass=task_class,
                    task_args=[source],
                    debug = debug))
            
        return processors,task_class,pause_qus

    def __report_file_names__(self,output_paths):
        print "\n[Status] This run will output the following files:"
        for src,subset_paths in output_paths.items():
            for subset,output_path in subset_paths.items():
                if not subset == "index":
                    print "\t%s" % (output_path.split("/")[-1],)
        print ""

    def __move_output_files__(self,output_paths):
        final_files = []
        for src, subset_paths in output_paths.items():
            for subset,output_path in subset_paths.items():
                if not subset == "index":
                    try:
                        move_location = output_path.replace(self._temp_dir,".")
                        shutil.move(output_path,move_location) #./ being the current working dir
                        final_files.append(move_location) 
                    except shutil.Error,e:
                        alt_filnm = "./%s_%s_%d.bam" % (src,subset,time.time()) 
                        print "[Warning] Output file may already exist, you may not" \
                        "have correct permissions for this file"
                        print "[Update]Trying to create output using unique filename:"
                        print "\t\t%s" % (alt_filnm,)
                        shutil.move(output_path,alt_filnm)
                        final_files.append(alt_filnm)
        return final_files

    def get_parser(self):
        #argparse imported in ./interface/parabam
        parser = self.default_parser()

        parser.add_argument('--output','-o',metavar='OUTPUT', nargs='+',required=False
            ,help="The name of the output that we wish to create. Must be same amount of space"\
            " separated entries as INPUT.")
        parser.add_argument('--debug',action="store_true",default=False,
            help="Only the first 5million reads will be processed")
        parser.add_argument('--pair',action="store_true",default=False
            ,help="A pair processor is used instead of a conventional processor")
        parser.add_argument('-r','--region',type=str,metavar="REGION",nargs='?',default=None
            ,help="The subset process will be run only on reads from this region\n"\
            "Regions should be colon seperated as specified by samtools (eg \'chr1:1000,5000\')")
        parser.add_argument('-s',type=int,metavar="INT",nargs='?',default=2
            ,help="Further parralise subset by running this many samples side-by-side. [Default 2]")
        parser.add_argument('-d',action="store_true",default=False,
            help="parabam will process reads marked PCR duplicate. Including this parameter may crash pair processing")
        parser.add_argument('-u',action="store_true",default=False,
            help="The files will be appended with a unique identifier. In the format <input_name>_<unique_id>_<subset_type>.bam")
        parser.add_argument('--counts',action="store_true",default=False,
            help="The amount of reads in each subset will be output as a .CSV alongside the .BAM file.")
        parser.add_argument('-v', choices=[0,1,2],default=0,type=int,
            help="Indicate the amount of information output by the program:\n"\
            "\t0: No output [Default]\n"\
            "\t1: Total Reads Processsed\n"\
            "\t2: Detailed output")

        return parser 

#...happily ever after
