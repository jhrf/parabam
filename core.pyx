import os
import time
import datetime
import sys
import Queue as Queue2
import gc
import shutil
import argparse

import pysam

import multiprocessing
from abc import ABCMeta, abstractmethod

#Tasks are started by parabam.Processor. Once started 
#they will carryout a predefined task on the provided
#subset of reads (task_set).
class Task(multiprocessing.Process):

    __metaclass__ = ABCMeta

    def __init__(self,object parent_bam,object task_set,
                     object outqu,int curproc,
                     object const,str parent_class):
        
        multiprocessing.Process.__init__(self)
        self._task_set = task_set
        self._outqu = outqu
        self._curproc = curproc
        self._temp_dir = const.temp_dir
        self._parent_class = parent_class
        self._parent_bam = parent_bam
        self.const = const

    def run(self):
        cdef dict results = self.__generate_results__()
        results["total"] = len(self._task_set)

        self._outqu.put(CorePackage(results=results,
                                destroy=self._destroy,
                                curproc=self._curproc,
                                parent_class=self._parent_class))
        #Trying to control mem useage
        del self._task_set
        gc.collect()

    def __get_temp_path__(self,typ,ext=".bam"):
        #self.pid ensures that the temp names are unique.
        return "%s/%s_%s_parabam_temp%s" % (self._temp_dir,typ,self.pid,ext)

    #Generate a dictionary of results using self._task_set
    @abstractmethod
    def __generate_results__(self):
        #Should always return a dictionary
        pass

cdef class Handler:

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object const,object pause_qu,dict out_qu_dict,object report=True):

        self._inqu = inqu
        self._pause_qu = pause_qu
        self._out_qu_dict = out_qu_dict

        self._report = report
        self._stats = {}

        self.const = const
        self._verbose = const.verbose
        self._output_paths = output_paths

        self._periodic_interval = 10

        self._destroy = False
        self._finished = False

        if const.verbose == 1:
            self._verbose = True
            self._update_output = self.__level_1_output__
        elif const.verbose == 2 or const.verbose == True:
            #In the case verbose is simply "True" or "level 2"
            self._verbose = True
            self._update_output = self.__level_2_output__
        else:#catching False and -v0
            self._verbose = False
            self._update_output = self.__level_2_output__

    def __standard_output__(self,out_str):
        sys.stdout.write(out_str + "\n")
        sys.stdout.flush()

    def __level_1_output__(self,out_str):
        #BUG! the fact this makes a call to __total_reads__ is ridiculous
        #this is making calls to a sub class method and just expecting it to be there.

        total_procd = self.__total_reads__()
        time = out_str.partition("Time: ")[2]
        sys.stdout.write("[Update] Processed: %d Time: %s\n" % (total_procd,time))
        sys.stdout.flush()

    def __level_2_output__(self,outstr):
        sys.stdout.write("\r" + outstr)
        sys.stdout.flush()

    #Must be overwritten if stats architecture is modififed
    def __total_reads__(self):
        return self._stats["total"]

    def listen(self,update_interval):
        destroy = self._destroy

        cdef int iterations = 0
        cdef int start_time = time.time()
        cdef int dealt  = 0
        cdef int curproc = 0

        update_output = self._update_output #speedup alias
        periodic_interval = self._periodic_interval #speedup alias
        finished = self.__is_finished__

        while not finished():
            iterations += 1
            #Listen for a process coming in...
            try:
                new_package = self._inqu.get(False)
                if new_package.destroy:
                    self._destroy = True

                if not new_package.results == {}:#If results are present...
                    self.__new_package_action__(new_package) #Handle the results
                    dealt += 1
                    if type(new_package) == CorePackage:
                        curproc = new_package.curproc

            except Queue2.Empty:
                #Queue empty. Continue with loop
                time.sleep(1)

            if iterations % periodic_interval == 0: 
                self.__periodic_action__(iterations)

            if self._verbose and self._report and iterations % update_interval == 0:
                outstr = self.__format_update__(curproc,start_time)
                update_output(outstr)

        self._inqu.close()
        self.__handler_exit__()

    def __is_finished__():
        if not self._destroy or not self._finished:
            return False
        return True

    def __format_update__(self,curproc,start_time):
        stats = []

        for stat_name in self._stats:
            if type(self._stats[stat_name]) is dict: #Account for divided stats.
                stat_update_str = "[%s] " % (stat_name,)
                for data in self._stats[stat_name]:
                    stat_update_str += "%s:%d " % (data,self._stats[stat_name][data])
                stat_update_str = stat_update_str[:-1] + " "

                stats.append(stat_update_str)
            else:
                stats.append("%s: %d" % (stat_name,self._stats[stat_name]))

        statstr = " ".join(stats)
        return "\r%s | Tasks: %d Time: %d " %\
                (statstr,curproc,self.__secs_since__(start_time),)

    #This code is a little ugly. Essentially, given a results
    #dictionary, it will go through and create a sensible output
    #string.
    def __auto_handle__(self,results,stats):
        for key_one in results:
            if type(results[key_one]) is int:
                if key_one in stats:
                    stats[key_one] += results[key_one]
                else:
                    stats[key_one] = results[key_one]
            elif type(results[key_one]) is dict:
                for key_two in results[key_one]:
                    if type(results[key_one][key_two]) is int:
                        if key_two in stats:
                            stats[key_two] += results[key_one][key_two]
                        else:
                            stats[key_two] = results[key_one][key_two]

    def __secs_since__(self,since):
        return int(time.time() - since)

    def __periodic_action__(self,iterations):
        #Overwrite with something that needs
        #to be done occasionally.       
        pass

    def __new_package_action__(self,new_package,**kwargs):
        #Handle the output of a task. Input will always be of type
        #parabam.Package. New results action should always call
        #self.__auto_handle__(new_package.results)
        pass

    def __handler_exit__(self,**kwargs):
        #When the handler finishes, do this
        pass

#The Processor iterates over a BAM, subsets reads and
#then starts a parbam.Task process on the subsetted reads
cdef class Processor:

    cdef int _chunk,_proc,_verbose

    def __init__(self,object parent_bam,object outqu,object const,
                 object TaskClass,list task_args,object pause_qu,
                 str region = None,object debug=False):

        self._debug = debug
        self.const = const
        self._chunk = const.chunk
        self._proc = const.proc
        self._verbose = const.verbose
        self._temp_dir = const.temp_dir

        self._region = region
        self._parent_bam = parent_bam
        self._parent_path = parent_bam.filename
        
        self._outqu = outqu

        #Class which we wish to run on each read
        self._TaskClass = TaskClass
        self._task_args = task_args 

        self._active_tasks = []
        self._active_count = 0
        self._prev = 0

    #Find data pertaining to assocd and all reads 
    #and divide pertaining to the chromosome that it is aligned to
    def run(self,int update):
    
        self.__pre_processor__()
        
        cdef object parent_object = self.__get_parent_object__(self._parent_path)
        cdef list collection = []
        cdef int collection_count = 0

        #function alias, for optimisation
        add_to_collection = self.__add_to_collection__
        start_task = self.__start_task__
        wait_for_tasks = self.__wait_for_tasks__

        #Insert a debug iterator into the processor, this workaround doesn't
        #slow processing when not in debug mode.
        parent_generator = self.__get_next_alig__(parent_object) if not self._debug \
                            else self.__get_next_alig_debug__(parent_object)

        for i,alig in enumerate(parent_generator):

            add_to_collection(parent_object,alig,collection)
            collection_count += 1

            if collection_count == self._chunk:
                wait_for_tasks(self._active_tasks,self._proc) # -2 for proc and handler 
                start_task(collection)      
                del collection
                collection = [] # already running
                collection_count = 0

        wait_for_tasks(self._active_tasks,0)
        start_task(collection,destroy=True)
        self.__end_processing__(parent_object)

    #__wait_for_tasks__ controls the amount of currently running processors
    #it waits until a processor is free and then allows the creation of a new
    #process
    def __wait_for_tasks__(self,list active_tasks,int max_tasks):
        update_tasks = self.__update_tasks__ #optimising alias
        update_tasks(active_tasks)
        cdef int currently_active = len(active_tasks)
        self._active_count = currently_active

        if max_tasks > currently_active:
            return

        while(max_tasks <= currently_active and currently_active > 0):
            update_tasks(active_tasks)
            currently_active = len(active_tasks)
            time.sleep(1)

        return

    def __update_tasks__(self,active_tasks):
        terminated_procs = []
        for pr in active_tasks:
            if not pr.is_alive():
                pr.terminate()
                terminated_procs.append(pr)
        
        for pr in terminated_procs:
            active_tasks.remove(pr)

        if len(terminated_procs) > 0:
            #Only invoked the GC if there is the possibility of
            #collecting any memory.
            del terminated_procs
            gc.collect()
                   
    def __start_task__(self,collection,destroy=False):
        pause = self.__query_pause_qu__(False)
        if pause:
            while True:
                pause = self.__query_pause_qu__(True)
                if not pause:
                    break
        
        args = {"task_set":collection,
                "outqu":self._outqu,
                "curproc":self._active_count,
                "parent_bam":self._parent_bam,
                "const":self.const,
                "parent_class"self.__class__.__name__}
        args.update(self._task_args)

        task = self._TaskClass(**args)
        task.start()
        self._active_tasks.append(task)

    def __get_master_bam__(self,master_file_path):
        if self.const.input_is_sam:
            return pysam.Samfile(self.const.master_file_path[self._source],"r")
        else:
            return pysam.Samfile(self.const.master_file_path[self._source],"rb")

    #Should be over written if we don't want data from BAM file.
    def __get_next_alig__(self,master_bam):
        if not self._region:
            for alig in master_bam.fetch(until_eof=True):
                yield alig
        else:
            for alig in master_bam.fetch(region=self._region):
                yield alig

    def __get_next_alig_debug__(self,master_bam):
        for i,alig in enumerate(master_bam.fetch(until_eof=True)):
            if i < 10000000:
                yield alig
            else:
                return

    def __end_processing__(self,master):
        master.close()
        for proc in self._active_tasks:
            proc.join()
            proc.terminate()

    def __pre_processor__(self):
        pass

    def __add_to_collection__(self,master,item,collection):
        collection.append(item)

    def __get_parent_bam__(self, master_file_path):
        return ParentAlignmentFile(master_file_path,self.const.input_is_sam)

     def __query_pause_qu__(self,bypass):
        #This fairly obfuscated code is to fix
        #a race condition when a pause can
        #be missed due to waiting on a previous pause
        try:
            if bypass:
                pause = True
            else:
                pause = self._inqu.get(False)

            if pause:
                while True:
                    try:
                        pause = self._inqu.get(False)
                        if not pause:
                            break
                    except Queue2.Empty:
                        time.sleep(1)
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

class Leviathon(object):
    #Leviathon takes objects of processors and handlers and
    #chains them together.
    def __init__(self,int max_processors, object const,list 
                processor_bundle,list handler_bundle,list queue_names,
                int update=10000000):

        self._max_processors = max_processors
        self._sources = sources
        self._processors = processors 
        self._handlers  = handlers
        self._update = update
        self._queue_names = queue_names
    
    def run(self,input_path,output_paths):
        parent = ParentAlignmentFile(input_path)
        references = parent.references
        lengths = parent.lengths

        chrm_by_length = reversed(sorted(zip(references,lengths),key=lambda tup: tup[1]))
        active_processors = []
        proccesor_pause_qus = []

        max_processors = self.max_processors

        default_qus = self.__create_queues__(self._queue_names) 

        handlers = self.__create_handlers__(handler_bundle,default_qus)
        handler_processes = self.__start_handlers__(handlers)

        while True:

            self.__update_processors__(active_processors,proccesor_pause_qus)
            self.__pause_routine__(default_qus["pause"],proccesor_pause_qus)

            if len(active_processors) < max_processors:
                try:
                    chrm,length chrm_by_length.next()
                    proc_class = processor_bundle["class"]
                    del processor_bundle["class"]

                    new_pause = Queue()
                    processor_bundle["out_qu"] = default_qus["main"]
                    processor_bundle["pause_qu"] = new_pause
                    processor_bundle["parent_bam"] = parent
                    processor_bundle["region"] = "%s" % (chrm,)

                    proc = proc_class(**processor_bundle)

                    ppr = multiprocessing.Process(target=proc.run,args=(self._update,))
                    ppr.start()

                    proccesor_pause_qus.append(new_pause)
                    active_processors.append(ppr)

                except StopIteration:
                    break    

        self.__wait_and_destroy_processors__(active_processors,proccesor_pause_qus)
        self.__destory__handlers__(handler_processes,handler_inqus)

        del active_processors
        del proccesor_pause_qus
        del handler_processes
        del handler_inqus

        gc.collect()

    def __create_queues__(self,queue_names):
        queues = {}
        for name in queue_names:
            queues[name] = Queue()
        return queues

    def __wait_and_destroy_processors__(self,active_processors,proccesor_pause_qus):
        for p,pause_qu in izip(active_processors,proccesor_pause_qus):
            p.join()
            pause_qu.close()

    def __pause_routine__(self,in_pause_qu,proccesor_pause_qus):
        pause = self.__pause_query__(in_pause_qu,proccesor_pause_qus)
        if pause:
            while True:
                pause = self.__pause_query__(in_pause_qu,proccesor_pause_qus,bypass=True)
                if not pause:
                    break

    #TODO: similar to code in the Processor class
    def __pause_query__(self,in_pause_qu,proccesor_pause_qus,bypass=True):
        try:
            if bypass:
                pause = True
            else:
                pause = in_pause_qu.get(False)
                if pause:
                    for qu in proccesor_pause_qus:
                        qu.put(True)

            if pause:
                while True:
                    try:
                        pause = in_pause_qu.get(False)
                        if not pause:
                            for qu in proccesor_pause_qus:
                                qu.put(False)
            
                    except Queue2.Empty:
                        time.sleep(1)
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

    #TODO: modified from Processor
     def __update_processors__(self,active_tasks,proccesor_pause_qus):
        terminated_procs = []
        terminated_qus = []
        for pr,pause_qu in active_tasks,proccesor_pause_qus:
            if not pr.is_alive():
                pr.terminate()
                terminated_procs.append(pr)
                terminated_qus.append(pause_qu)

        for pr,qu in (terminated_procs,proccesor_pause_qus):
            active_tasks.remove(pr)
            
            qu.close()
            proccesor_pause_qus.remove(qu)

        if len(terminated_procs) > 0:
            #Only invoked the GC if there is the possibility of
            #collecting any memory.
            del terminated_procs
            del terminated_qus
            gc.collect()

    def __create_handlers__(self,handler_order,handler_bundle,queues,parent_bam,output_paths):
        handlers = []
        handler_inqus = []

        for handler_class in handler_order:
            handler_args = dict(hander_bundle[handler_class])

            handler_args["parent_bam"] = parent_bam
            handler_args["output_paths"] = output_paths

            #replace placeholder with queues
            handler_args["pause_qu"] = queues["pause"]
            handler_args["inqu"] = queues[handler_args["inqu"]]
            handler_args["out_qu_dict"] = dict(\
                    [(name,queues[name]) for name in handler_args["out_qu_dict"] ])

            handler_inqus.append(handler_args["inqu"])
            handlers.append(handler_class(**handler_args))

        return handlers,handler_inqus

    def __start_handlers__(self,handlers):
        handler_processes = []
        for handler in handlers:
            hpr = multiprocessing.Process(target=h.listen,args=(self._update,))
            hpr.start()
            handler_processes.append(hpr)
        return handler_processes

    def __destory__handlers__(self,handler_processes,handler_inqus):
        for process,queue in izip(handler_processes,handler_inqus):
            queue.put(DestroyPackage())
            process.join()

#Provides a conveinant way for providing an Interface to parabam
#programs. Includes default command_args and framework for 
#command-line and programatic invocation. 
class Interface(object):

    __metaclass__ = ABCMeta

    def __init__(self,temp_dir):
        self._temp_dir = temp_dir

    def __introduce__(self,name):
        intro =  "%s has started. Start Time: " % (name,)\
            + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        underline = "-" * len(intro)
        print intro
        print underline

    def __goodbye__(self,name):
        print "%s has finished. End Time: " % (name,)\
            + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

    def __sort_and_index__(self,fnm,verbose=False,tempDir=".",name=False):
        if not os.path.exists(fnm+".bai"):
            if verbose: print "%s is not indexed. Sorting and creating index file." % (fnm,)
            tempsort_path = self.__get_unique_temp_path__("SORT",temp_dir=self._temp_dir)
            optStr = "-nf" if name else "-f"
            pysam.sort(optStr,fnm,tempsort_path)
            os.remove(fnm)
            os.rename(tempsort_path,fnm)
            if not name: pysam.index(fnm)

    def __get_unique_temp_path__(self,temp_type,temp_dir="."):
        #TODO: Duplicated from parabam.interface.merger. Find a way to reformat code
        #to remove this duplication
        return "%s/%sTEMP%d.bam" % (temp_dir,temp_type,int(time.time()),)

    def __get_group__(self,bams,names=[],multi=1):
        if multi == 0:
            multi = 1
        if names == []:
            names = [self.__get_basename__(path) for path in bams]
        for i in xrange(0,len(bams),multi):
            yield (bams[i:i+multi],names[i:i+multi])

    def __get_basename__(self,path):
        base,ext = os.path.splitext(os.path.basename(path))
        return base

    def default_parser(self):

        parser = ParabamParser(conflict_handler='resolve',
                    formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-p',type=int,nargs='?',default=4
            ,help="The maximum amount of tasks run concurrently. This number should\n"\
            "be less than or equal to the amount of cores on your machine [Default: 4]")
        parser.add_argument('-c',type=int,nargs='?',default=100000
            ,help="How many reads each parallel task is given. A higher number\n"\
            "will mean faster analysis but also a greater burden on memory [Default:100000]")
        return parser

    @abstractmethod
    def run_cmd(self,parser):
        #This is usualy just a function that
        #takes an argparse parser and turns 
        #passes the functions to the run function
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def get_parser(self):
        pass

class ParabamParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help()
        sys.stderr.write('\nerror: %s\n' % message)
        sys.exit(2)

class Const(object):
    
    def __init__(self,temp_dir,verbose,chunk,proc,**kwargs):
        self.temp_dir = temp_dir
        self.verbose = verbose
        self.chunk = chunk
        self.proc = proc

        for key, val in kwargs.items():
            setattr(self,key,val)

    def add(self,key,val):
        setattr(self,key,val)

class Package(object):
    def __init__(self,results):
        self.results = results

class CorePackage(Package):
    def __init__(self,results,curproc,parent_class):
        super(CorePackage,self).__init__(results)
        self.curproc = curproc
        self.parent_class = parent_class

class DestroyPackage(Package):
    def __init__(self):
        super(DestroyPackage,self).__init__(results={})
        self.destroy = True

class ParentAlignmentFile(object):
    
    def __init__(self,path,input_is_sam=False):
        has_index = os.path.exists(os.path.join("%s%s" % (path,".bai")))

        if input_is_sam:
            mode = "r"
        else:
            mode = "rb"
        parent = pysam.AlignmentFile(path,mode)
        self.filename = parent.filename
        self.references = parent.references
        self.header = parent.header
        self.lengths = parent.lengths

        if not input_is_sam and has_index:
            self.mapped = parent.mapped
            self.nocoordinate = parent.nocoordinate
            self.nreferences = parent.nreferences
            self.unmapped = parent.unmapped
        else:
            self.mapped = 0
            self.nocoordinate = 0
            self.nreferences = 0
            self.unmapped = 0

        parent.close()

    def getrname(self,tid):
        return self.references[tid]

    def gettid(self,reference):
        for i,ref in enumerate(self.references):
            if reference == ref:
                return i
        return -1


#And they all lived happily ever after...
