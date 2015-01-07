import pysam
import pdb
import parabam
import time
import sys
import os
import gc
import shutil
import gzip

from itertools import izip
from collections import deque, Counter
from multiprocessing import Queue,Process
from abc import ABCMeta, abstractmethod
from parabam.support import HandlerMerge,MergePackage

class TaskSubset(parabam.core.Task):

    def __init__(self,object task_set,object outqu,object curproc,object destroy,object const,str source):
        super(TaskSubset, self).__init__(task_set=task_set,
                                        outqu=outqu,
                                        curproc=curproc*len(const.sources),
                                        destroy=destroy,
                                        const=const)
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
            else:
                sys.stdout.write("[ERROR] Unrecognised return type from user engine!\n")
                print subset_decision
                print type(subset_decision)
                print self.__class__
                sys.stdout.flush()

    def __engine__(self,read,user_constants,master):
        self.const.user_engine(read,user_constants,master)
        
class PairTaskSubset(TaskSubset):

    def __init__(self,object task_set,object outqu,object curproc,object destroy,object const,str source):
        super(PairTaskSubset, self).__init__(task_set=task_set,
                                        outqu=outqu,
                                        curproc=curproc*len(const.sources),
                                        destroy=destroy,
                                        const=const,
                                        source=source)

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

    def __read_pair_decision__(self,read1,read2,engine,user_constants,master):
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
        try:
            mate = loners[read.qname]
            del loners[read.qname]
            return read,mate

        except KeyError:
            loners[read.qname] = read
            return None,None

class HandlerSubset(parabam.core.Handler):

    def __init__(self,object inqu,object outqu,object const,object destroy_limit):
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

    def __get_total_processed_reads__(self):
        total = 0
        for source,deep_stat in self._stats.items():
            for name,value in deep_stat.items():
                if name == "total":
                    total += value
        return total

    def __new_package_action__(self,new_package,**kwargs):
        results = new_package.results
        source = results["source"]
        self.__auto_handle__(results,source)

        for subset in self._subset_types:
            self._merge_stores[source][subset].append((results["counts"][subset],results["temp_paths"][subset],))

    def __periodic_action__(self,iterations):
        for source in self._sources:
            for subset in self._subset_types:
                if self.__test_merge_store__(source,subset):
                    self.__add_merge_task__(name=self._output_paths[source][subset],
                                            results=self._merge_stores[source][subset],
                                            subset_type=subset,source=source,
                                            total=self._stats[source]["total"])
                    self.__update_merge_store__(source,subset)
                    
    def __update_merge_store__(self,source,subset):
        self._mergecount += 1
        self._merge_stores[source][subset] = []
        gc.collect()

    def __test_merge_store__(self,source,subset):
        return self._merge_stores[source][subset] > 0
            
    def __add_merge_task__(self,name,results,subset_type,source,total,destroy=False):
        
        res = MergePackage(name=name,results=list(results),
                            subset_type=subset_type,source=source,
                            destroy=destroy,total=total,time_added=time.time())
        self._mergequeue.put(res)

    def __total_reads__(self):
        return sum(map(lambda s : self._stats[s]["total"],self._sources))

    def __handler_exit__(self,**kwargs):
        if self._verbose:
            self.__standard_output__("\n[Status] Processed %d reads from bam files\n" % (self.__total_reads__(),))
            self.__standard_output__("[Status] Waiting for merge operation to finish...\n")

        for source in self._sources:
            for subset in self._subset_types:
                self.__add_merge_task__(name=self._output_paths[source][subset],
                                results=self._merge_stores[source][subset],subset_type=subset,
                                source=source,total=self._stats[source]["total"],
                                destroy=True)
                self.__update_merge_store__(source,subset)

class HandlerIndex(parabam.core.Handler):

    def __init__(self,object inqu,object mainqu,object const,object destroy_limit, object TaskClass):
        super(HandlerIndex,self).__init__(inqu,const,destroy_limit=destroy_limit,report=False)

        self._sources = const.sources
        self._subset_types = const.subset_types
        
        self._mainqu = mainqu
        self._TaskClass = TaskClass
        
        self._index_paths = self.__instalise_index_queues__()

        self._periodic_interval = 1

        self._total_loners = 0
        self._rescued = 0

    def __instalise_index_queues__(self):
        loner_types = ["both_mapped","both_unmapped","one_mapped"]
        index_queues = {}
        for source in self.const.sources:
            index_queues[source] = {}
            for loner_type in loner_types:
                index_queues[source][loner_type] = deque()
        return index_queues

    def __new_package_action__(self,new_package,**kwargs):
        subset_type = new_package.subset_type

        if subset_type == "index":
            self.__start_primary_task__(new_package)
        elif subset_type == "primary":
            self.__handle_primary_task__(new_package)
        elif subset_type == "match_maker":
            self.__handle_match_maker__(new_package)
 
    def __handle_match_maker__(self,new_package):
        source = new_package.source
        name = new_package.name

        for path in new_package.results:
            self._rescued += (new_package.total * 2) #this number is by pairs, to get read num we times 2
            self._index_paths[source][name].appendleft(path)

    def __handle_primary_task__(self,new_package):
        source = new_package.source
        for loner_type,loner_path in new_package.results.items():
            self._index_paths[source][loner_type].appendleft(loner_path)

    def __start_primary_task__(self,new_package):
        for loner_count,path in new_package.results:
            self._total_loners += loner_count
            new_task = PrimaryTask(path,self._inqu,self.const,new_package.source)
            new_task.start()

    def __periodic_action__(self,iterations):
        sys.stdout.write("\rRescued: %d Total Loners: %d" % (self._rescued,self._total_loners))
        sys.stdout.flush()
        for source,queue_dict in self._index_paths.items():
            for loner_type,queue in queue_dict.items():
                path_1 = None
                try:
                    for x in xrange(3):
                        path_1 = queue.pop()
                        path_2 = queue.pop()
                        
                        new_task = MatchMakerTask([path_1,path_2],self._inqu,self.const,source,loner_type,
                                            {"queue":self._mainqu,"const":self.const,
                                             "task_args":[source],"TaskClass":self._TaskClass})
                        new_task.start()
                        path_1 = None

                except IndexError:
                    if not path_1 == None:
                        queue.appendleft(path_1)
                        #If we fail on path1 we need to put it back.

    def __handler_exit__(self,**kwargs):
        print "Exiting..."

class ChaserClass(Process):
    def __init__(self):
        super(ChaserClass,self).__init__()

    def __sort_and_index__(self,path):
        temp_path = "%s/sort_temp_%d_%s" % (self.const.temp_dir,time.time(),self.pid)
        pysam.sort(path,temp_path)
        temp_path += ".bam"

        os.remove(path)
        os.rename(temp_path,path)

        pysam.index(path)

    def __read_generator__(self,alig_file_obj):
        for read in alig_file_obj.fetch(until_eof=True):
            yield read

    def __get_loner_temp__(self,loner_type,uid=0):
        return "%s/%s_%s_%d_%d.bam" % (self.const.temp_dir,loner_type,self.pid,time.time(),uid,)


class PrimaryTask(ChaserClass):

    def __init__(self,object unsorted_path,object outqu,object const,str source):
        super(PrimaryTask,self).__init__()
        
        self._unsorted_path = unsorted_path
        self._outqu = outqu
        self._source = source
        self.const = const

    def run(self):
        #speedup
        write_loner = self.__write_loner_type__
        #speedup

        unsorted_path = self._unsorted_path

        self.__sort_and_index__(unsorted_path)
        unsorted_object = pysam.AlignmentFile(unsorted_path,"rb")

        loner_types = ["both_mapped","both_unmapped","one_mapped"]
        loner_paths = self.__create_loner_paths__(loner_types)
        loner_objects = self.__create_loner_objects__(loner_paths,unsorted_object)

        loner_count = 0

        for read in self.__read_generator__(unsorted_object):
            current_type = ""
            if read.is_unmapped and read.mate_is_unmapped:
                current_type = "both_unmapped"
            elif not read.is_unmapped and not read.mate_is_unmapped:
                current_type = "both_mapped"
            else:
                current_type = "one_mapped"

            write_loner(loner_objects,current_type,read)
            loner_count += 1

        self.__close_loner_objects__(loner_objects)
        os.remove(unsorted_path)
        os.remove(unsorted_path+".bai")

        loner_pack = MergePackage(name="",results=loner_paths,
                subset_type="primary",source=self._source,
                destroy=False,total=loner_count,time_added=time.time())
        self._outqu.put(loner_pack)

    def __create_loner_objects__(self,loner_paths,unsorted_object):
        objects = {}
        for loner_type,path in loner_paths.items():
            objects[loner_type] = pysam.AlignmentFile(path,"wb",template=unsorted_object)
        return objects

    def __create_loner_paths__(self,loner_types):
        paths = {}
        for loner_type in loner_types:
            paths[loner_type] = self.__get_loner_temp__(loner_type)
        return paths

    def __write_loner_type__(self,loner_objects,loner_type,read):
        loner_objects[loner_type].write(read)

    def __close_loner_objects__(self,loner_objects):
        for loner_type,loner_object in loner_objects.items():
            loner_object.close()

class MatchMakerTask(ChaserClass):

    def __init__(self,object index_paths,object outqu,object const,str source,str loner_type,object child_package):
        super(MatchMakerTask,self).__init__()
        
        self._child_package = child_package
        self._index_paths = index_paths
        self._outqu = outqu
        self._source = source
        self._loner_type = loner_type
        self.const = const

    def run(self):
        path1,path2 = self._index_paths
        loner_path = self.__get_loner_temp__(self._loner_type)

        path1_object = self.__create_path_object__(path1)
        path2_object = self.__create_path_object__(path2)
        loner_object = pysam.AlignmentFile(loner_path,"wb",template=path1_object)

        pairs,loners = self.__instalise_loners_and_pairs__(path1_object)

        for read in self.__read_generator__(path2_object):
            read1,read2 = self.__query_loners__(read,loners)
            if read1:
                pairs[read1.qname] = (read1,read2)

        stash_count,loner_paths = self.__stash_loners__(loners,self.const.chunk)
        map( lambda fil: fil.close() , [path1_object,path2_object])

        if len(loners) > 0:
            results = [loner_path]
        else:
            results = []
            os.remove(loner_path)

        if len(pairs) > 0:
            #Send pairs onto the subset decision
            self.__launch_child_task__(pairs)

        loner_pack = MergePackage(name=self._loner_type,results=results,
                    subset_type="match_maker",source=self._source,
                    destroy=False,total=len(pairs),time_added=time.time())
        self._outqu.put(loner_pack)

        del loners        
        os.remove(path1)
        os.remove(path2)
        os.remove(path1+".bai")
        os.remove(path2+".bai")
        gc.collect()

    def __create_path_object__(self,path):
        self.__sort_and_index__(path)
        path_object = pysam.AlignmentFile(path,"rb")
        return path_object

    def __launch_child_task__(self,pairs):
        child_pack = self._child_package
        args = [pairs,child_pack["queue"],0,False,child_pack["const"]]
        args.extend(child_pack["task_args"])
        task = child_pack["TaskClass"](*args)
        task.start()

    def __stash_loners__(self,loners,loner_object,chunk,template_object):
        loner_count = len(loners)
        path_count = (loner_count // chunk) + 1
        paths = self.__get_loner_paths__(path_count)
        loner_objects = self.__get_loner_objects__()
        for i,(qname,read) in enumerate(loners.items()):
            file_num = i % 3
            loner_objects[file_num].write(read)

        map(lambda filnm : fil.close() , loner_objects) #close files
        return stashed,paths

    def __get_loner_paths__(self,amount):
        paths = []
        for n in xrange(amount):
            paths.append(self.__get_loner_temp__(self._loner_type,n))
        return paths

    def __get_loner_objects__(self,paths,template_object):
        objects = []
        for path in paths:
            pysam.AlignmentFile(path,"wb",template=template_object)
        return objects

    def __query_loners__(self,read,loners):
        try:
            mate = loners[read.qname]
            del loners[read.qname]
            return read,mate
        except KeyError:
            loners[read.qname] = read
            return None,None

    def __instalise_loners_and_pairs__(self,alig_file_obj):
        query_loners = self.__query_loners__
        loners = {}
        pairs = {}
        for read in alig_file_obj.fetch():
            read1,read2 = self.__query_loners__(read,loners)
            if read1:
                pairs[read1.qname] = (read1,read2) 
        return pairs,loners

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

class PairProcessor(ProcessorSubset):

    def __init__(self,object outqu,object const,object TaskClass,object task_args,object debug=False):
        super(PairProcessor,self).__init__(outqu,const,TaskClass,task_args,debug)
        
    def __add_to_collection__(self,master,item,collection):
        collection.append(item)

class Interface(parabam.core.UserInterface):
    """The interface to parabam subset.
    Users will primarily make use of the ``run`` function."""

    def __init__(self,temp_dir,exe_dir):
        super(Interface,self).__init__(temp_dir,exe_dir)

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
            outputs= cmd_args.output,
            proc= cmd_args.p,
            chunk= cmd_args.c,
            verbose= verbose,
            subset_types= subset_types,
            user_constants = user_constants,
            user_engine = user_engine,
            engine_is_class = False,
            fetch_region = cmd_args.region,
            pair_process=cmd_args.m,
            side_by_side = cmd_args.s,
            debug = cmd_args.debug
            )
    
    def run(self,input_bams,outputs,proc,chunk,subset_types,
            user_constants,user_engine,fetch_region=None,side_by_side=2,
            keep_in_temp=False,engine_is_class=False,verbose=False,
            pair_process=False,debug=False):

        if not outputs or not len(outputs) == len(input_bams):
            print "[Status] Using default naming scheme."
            outputs = [ self.__get_basename__(b) for b in input_bams ]

        #AT SOME POINT WE SHOULD HANDLE UNSORTED BAMS. EITHER HERE OR AT THE PROCESSOR
        final_files = []

        for input_group,output_group in self.__get_group__(input_bams,outputs,multi=side_by_side):
            
            output_paths = dict([(source,{}) for source in output_group])
            master_file_path = {}

            if pair_process:
                #This index subset records the failed attempts to find mates
                #if we are pair processing
                subset_types.append("index")

            for mst,source in zip(input_group,output_group):
                master_file_path[source] = mst
                for typ in subset_types:
                    output_paths[source][typ] = "%s/%s_%s.bam" % (self._temp_dir,source.replace(".bam",""),typ,)
                    
            if verbose: self.__report_file_names__(output_paths)

            const = parabam.core.Const(output_paths=output_paths,
                                temp_dir=self._temp_dir,
                                master_file_path=master_file_path,
                                chunk=chunk,proc=(proc // len(input_group)),
                                verbose=verbose,thresh=0,
                                subset_types=subset_types,
                                sources=output_group,
                                exe_dir=self._exe_dir,
                                user_constants=user_constants,
                                user_engine=user_engine,
                                fetch_region=fetch_region,
                                pair_process=pair_process)

            task_qu = Queue()
            processors,task_class = self.__create_processors__(task_qu,const,debug,engine_is_class)
            handlers = self.__create_handlers__(task_qu,const,task_class)

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
                            final_files.append(path)
            else:
                final_files.extend(self.__move_output_files__(output_paths))
            
            gc.collect()

        return final_files

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

    def __create_handlers__(self,task_qu,object const,task_class):
        handlers = []
        merge_qu = Queue()
        merge_out = Queue()

        handlers.append(HandlerSubset(inqu=task_qu,outqu=merge_qu,
                        const=const,destroy_limit=len(const.sources)))

        handlers.append(HandlerMerge(inqu=merge_qu,outqu=merge_out,
                        const=const,destroy_limit=len(const.sources)*len(const.subset_types)))

        if const.pair_process:
            handlers.append(HandlerIndex(inqu=merge_out,mainqu=task_qu,
                        const=const,destroy_limit=len(const.sources),TaskClass=task_class))

        return handlers

    def __create_processors__(self,task_qu,object const,debug,engine_is_class):
        processors = []

        if const.pair_process:
            processor_class = PairProcessor
            task_class = PairTaskSubset
        else:
            processor_class = ProcessorSubset
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
            processors.append(processor_class(outqu=task_qu,
                const=const,
                TaskClass=task_class,
                task_args=[source],
                debug = debug))
        
        return processors,task_class

    def __report_file_names__(self,output_paths):
        print "[Status] This run will output the following files:"
        for src,subset_paths in output_paths.items():
            for subset,output_path in subset_paths.items():
                print "\t%s" % (output_path.split("/")[-1],)
        print ""

    def __move_output_files__(self,output_paths):
        final_files = []
        for src, subset_paths in output_paths.items():
            for subset,output_path in subset_paths.items():
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

        parser.add_argument('-r','--region',type=str,metavar="REGION",nargs='?',default=None
            ,help="The subset process will be run only on reads from this region\n"\
            "Regions should be colon seperated as specified by samtools (eg \'chr1:1000,5000\')")
        parser.add_argument('--output','-o',metavar='OUTPUT', nargs='+',required=False
            ,help="The name of the output that we wish to create. Must be same amount of space"\
            " separated entries as INPUT.")
        parser.add_argument('-s',type=int,metavar="INT",nargs='?',default=2
            ,help="Further parralise subset by running this many samples side-by-side. [Default 2]")
        parser.add_argument('--debug',action="store_true",default=False,
            help="Only the first 5million reads will be processed")
        parser.add_argument('-m',action="store_true",default=False
            ,help="A pair processor is used instead of a conventional processor")
        parser.add_argument('-v', choices=[0,1,2],default=0,type=int,
            help="Indicate the amount of information output by the program:\n"\
            "\t0: No output [Default]\n"\
            "\t1: Total Reads Processsed\n"\
            "\t2: Detailed output")

        return parser 

#...happily ever after