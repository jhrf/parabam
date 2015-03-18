import os
import time
import datetime
import sys
import Queue as Queue2
import gc
import shutil
import argparse

import pysam

from itertools import izip
from multiprocessing import Process,Queue
from abc import ABCMeta, abstractmethod

#Tasks are started by parabam.Processor. Once started 
#they will carryout a predefined task on the provided
#subset of reads (task_set).
class TTask(Process):

    __metaclass__ = ABCMeta

    def __init__(self,object parent_bam,object task_set,
                     object outqu,int curproc,
                     object const,str parent_class):
        
        Process.__init__(self)
        self._task_set = task_set
        self._outqu = outqu
        self._curproc = curproc
        self._temp_dir = const.temp_dir
        self._parent_class = parent_class
        self._parent_bam = parent_bam
        self.const = const

    def run(self):
        #print "Task Started %d" % (self.pid,)
        results = self.__generate_results__()
        #print "Task Finished %d" % (self.pid,)
        results["total"] = len(self._task_set)

        self._outqu.put(CorePackage(results=results,
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

        self._parent_bam = parent_bam
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
        if self._stats == {}:
            return 0
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
                if type(new_package) == DestroyPackage:
                    self._destroy = True

                if not new_package.results == {}:#If results are present...
                    self.__new_package_action__(new_package) #Handle the results
                    dealt += 1
                    if type(new_package) == CorePackage:
                        curproc = new_package.curproc

            except Queue2.Empty:
                #Queue empty. Continue with loop
                time.sleep(.5)

            if iterations % periodic_interval == 0: 
                self.__periodic_action__(iterations)

            if self._verbose and self._report and iterations % update_interval == 0:
                outstr = self.__format_update__(curproc,start_time)
                update_output(outstr)

        self._inqu.close()
        self.__handler_exit__()

    def __is_finished__(self):
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

class Task(Process):

    def __init__(self,parent_bam,seek,mod,outqu,size,header,temp_dir):
        super(Task,self).__init__()
        self._parent_bam = parent_bam
        self._seek = seek
        self._mod = mod
        self._outqu = outqu
        self._size = size
        self._header = header
        self._temp_dir = temp_dir

    def run(self):
        cdef int size = self._size
        cdef int sub_count = 0

        temp = pysam.AlignmentFile("%s/%s_%d_%d_subset.bam" % (self._temp_dir,self._mod,self.pid,time.time()),"wb",header = self._header)
        bamfile = pysam.AlignmentFile(self._parent_bam.filename,"rb")
        bamiter = bamfile.fetch(until_eof=True)
        bamfile.seek(self._seek)

        next_read = bamiter.next

        for i in xrange(size):
            try:
                read = next_read()
                if "TTAGGGTTAGGG" in read.seq:
                    temp.write(read)
                    sub_count += 1
            except StopIteration:
                break

        results = {"temp_paths":{"subset":temp.filename},
                   "counts": {"subset":sub_count},
                   "total" : size}
        temp.close()
        self._outqu.put(CorePackage(results=results,
                        curproc=self._mod,
                        parent_class=self.__class__.__name__))

#The Processor iterates over a BAM, subsets reads and
#then starts a parbam.Task process on the subsetted reads
class Processor(Process):
    def __init__(self,str input_path,int mod,object outqu,object const,int size):
        super(Processor,self).__init__()

        self._mod = mod
        self._input_path = input_path
        self._proc = const.proc
        self._temp_dir = const.temp_dir
        self._outqu = outqu
        self._size = size

    #Find data pertaining to assocd and all reads 
    #and divide pertaining to the chromosome that it is aligned to
    def run(self):

        #Insert a debug iterator into the processor, this workaround doesn't
        #slow processing when not in debug mode.

        parent_bam = pysam.AlignmentFile(self._input_path,"rb")
        parent_iter = parent_bam.fetch(until_eof=True)
        parent_generator = self.__bam_generator__(parent_iter)
        
        for i,command in enumerate(parent_generator):
            task = Task(parent_bam,parent_bam.tell(),self._mod,self._outqu,1000000,parent_bam.header,self._temp_dir)
            task.start()
        parent_bam.close()
        print "generator finished %d" % (self._mod,)

    def __bam_generator__(self,parent_iter):
        cdef int mod = self._mod
        cdef int proc = self._proc
        cdef int iterations = 0
        cdef int chunk = 1000000

        while True:            
            try:
                if iterations % proc == mod:
                    yield True
                    for i in xrange(chunk):
                        parent_iter.next()
                else:
                    for i in xrange(chunk):
                        parent_iter.next()
                iterations += 1
            except StopIteration:
                break
        return
        
class Leviathon(object):
    #Leviathon takes objects of processors and handlers and
    #chains them together.
    def __init__(self,int max_processors, object const,dict processor_bundle,
                dict handler_bundle,list handler_order,list queue_names,int update,str region):

        self._max_processors = max_processors
        self._processor_bundle = processor_bundle 
        self._handler_bundle  = handler_bundle
        self._handler_order = handler_order
        self._region = region
        self._update = update
        self._queue_names = queue_names
    
    def run(self,input_path,output_paths):
        parent = ParentAlignmentFile(input_path)

        default_qus = self.__create_queues__(self._queue_names) 

        processor_bundle = dict(self._processor_bundle)
        self.__run_info_to_bundle__(processor_bundle,default_qus,parent)

        handlers,handler_inqus = self.__create_handlers__(self._handler_order,self._handler_bundle,
                                                          default_qus,parent,output_paths)
        handler_processes = self.__start_handlers__(handlers)

        max_processors = self._max_processors
        proccesor_pause_qus = []
        active_processors = []

            #self.__update_processors__(active_processors,proccesor_pause_qus)
            #self.__pause_routine__(default_qus["pause"],proccesor_pause_qus)

        for submit_size,proc_num in izip(self.__power_generator__(),self.__proc_num_generator__(self._max_processors)):
            proc_proc = Processor(input_path,proc_num,default_qus["main"],self._processor_bundle["const"],2*(10**submit_size))
            proc_proc.start()
            active_processors.append(proc_proc)

        for proc_proc in active_processors:
            proc_proc.join()

        self.__destory__handlers__(handler_processes,handler_inqus)

        for hand_proc in handler_processes:
            hand_proc.join()

        del active_processors
        del proccesor_pause_qus
        del handler_inqus

        gc.collect()

    def __proc_num_generator__(self,total_proc):
        for i in reversed(xrange(total_proc)):
            yield i

    def __power_generator__(self):
        while True:
            yield 5
            yield 6

    def __run_info_to_bundle__(self,processor_bundle,default_qus,parent):
        processor_bundle["outqu"] = default_qus[processor_bundle["outqu"]]
        processor_bundle["parent_bam"] = parent.filename

    def __create_queues__(self,queue_names):
        queues = {}
        for name in queue_names:
            queues[name] = Queue()
        queues["pause"] = Queue()
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
    def __pause_query__(self,in_pause_qu,proccesor_pause_qus,bypass=False):
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
                            break
                    except Queue2.Empty:
                        time.sleep(.5)
        except Queue2.Empty:
            pass

        last = False  
        while True:
            try:
                pause = in_pause_qu.get(False)
                last = pause
            except Queue2.Empty:
                break
        return last

    #TODO: modified from Processor
    def __update_processors__(self,active_tasks,proccesor_pause_qus):
        terminated_procs = []
        terminated_qus = []
        for pr,pause_qu in izip(active_tasks,proccesor_pause_qus):
            if not pr.is_alive():
                pr.terminate()
                terminated_procs.append(pr)
                terminated_qus.append(pause_qu)

        for pr,qu in izip(terminated_procs,terminated_qus):
            active_tasks.remove(pr)

            qu.close()
            proccesor_pause_qus.remove(qu)

        if len(terminated_procs) > 0:
            #Only invoked the GC if there is the possibility of
            #collecting any memory.
            del terminated_procs
            del terminated_qus
            gc.collect()

    def __create_handlers__(self,handler_order,handler_bundle,
                            queues,parent_bam,output_paths):
        handlers = []
        handler_inqus = []

        for handler_class in handler_order:
            handler_args = dict(handler_bundle[handler_class])

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
            hpr = Process(target=handler.listen,args=(self._update,))
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
