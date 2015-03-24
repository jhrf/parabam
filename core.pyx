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

cdef class Handler:

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object constants,object pause_qus,dict out_qu_dict,object report=True):

        self._parent_bam = parent_bam
        self._inqu = inqu
        self._pause_qus = pause_qus
        self._out_qu_dict = out_qu_dict

        self._report = report
        self._stats = {}

        self._constants = constants
        self._verbose = constants.verbose
        self._output_paths = output_paths

        self._periodic_interval = 10

        self._destroy = False
        self._finished = False

        if constants.verbose == 1:
            self._verbose = True
            self._update_output = self.__level_1_output__
        elif constants.verbose == 2 or constants.verbose == True:
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

            except Queue2.Empty:
                #Queue empty. Continue with loop
                time.sleep(.5)

            if iterations % periodic_interval == 0: 
                self.__periodic_action__(iterations)

            if self._verbose and self._report and iterations % update_interval == 0:
                outstr = self.__format_update__(start_time)
                update_output(outstr)

        self._inqu.close()
        self.__handler_exit__()

    def __is_finished__(self):
        if not self._destroy or not self._finished:
            return False
        return True

    def __format_update__(self,start_time):
        stats = "[Update] "
        for name,value in self._stats.items():
            stats += "%s:%d " % (name,value)
        return "\r%s| Time: %d " %\
                (stats,self.__secs_since__(start_time),)

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
    __metaclass__=ABCMeta

    def __init__(self,parent_bam,inqu,outqu,task_size,constants):
        super(Task,self).__init__()
        self._parent_bam = parent_bam
        self._parent_path = self._parent_bam.filename
        self._header = self._parent_bam.header

        self._outqu = outqu
        self._inqu = inqu

        self._task_size = task_size
        self._temp_dir = constants.temp_dir
        self._constants = constants 
        
        self._dealt = 0

    def run(self):
        bamfile = pysam.AlignmentFile(self._parent_path,"rb")
        iterator = bamfile.fetch(until_eof=True)
        next_read = iterator.next

        generate_results = self.__generate_results__

        while True:
            try:
                package = self._inqu.get(False)
                
                if type(package) == DestroyPackage:
                    bamfile.close()
                    del iterator
                    del bamfile
                    break

                seek = package                
                bamfile.seek(seek)
                time.sleep(.01)

                results = generate_results(iterator)
                
                self._dealt += 1
                time.sleep(0.005)

                self._outqu.put(Package(results=results))

            except Queue2.Empty:
                time.sleep(5)
            except StopIteration:
                results = self.__handle_stop_iteration__(seek,bamfile,iterator)
                self._dealt += 1 #Copy paste code here.
                time.sleep(0.005) #And here
                self.__post_run_routine__()
                self._outqu.put(Package(results=results))

        return

    def __generate_results__(self,iterator,**kwargs):
        self.__pre_run_routine__(iterator)
        self.__process_task_set__(iterator)
        results = self.__get_results__()
        results["total"] = self._task_size
        self.__post_run_routine__()
        return results

    def __handle_stop_iteration__(self,seek,bamfile,iterator):
        cdef int count = 0 
        bamfile.seek(seek)
        try: #Count the amount of reads until end of file
            for read in iterator:
                count += 1
        except StopIteration:
            pass

        bamfile.seek(seek)

        try:
            self.__process_task_set__(iterator)
        except StopIteration:
            pass

        results = self.__get_results__()
        results["total"] = count
        return results

    def __get_temp_object__(self,path):
        return pysam.AlignmentFile(path,"wb",header=self._parent_bam.header)

    def __get_temp_path__(self,identity):
        file_name = "%s_%d_%d_%s" %\
            (identity,self.pid,self._dealt, os.path.split(self._parent_bam.filename)[1])
        return os.path.join(self._temp_dir,file_name)

    @abstractmethod
    def __pre_run_routine__(self,iterator,**kwargs):
        pass

    @abstractmethod
    def __process_task_set__(self,iterator,**kwargs):
        pass

    @abstractmethod
    def __get_results__(self):
        pass

    @abstractmethod
    def __post_run_routine__(self,**kwargs):
        pass

#The FileReader iterates over a BAM, subsets reads and
#then starts a parbam.Task process on the subsetted reads
class FileReader(Process):
    def __init__(self,str input_path,int proc_id,object outqu,
                 int task_n,object constants,object Task,
                 object pause_qu):
        
        #TODO FETCH REGION!
        super(FileReader,self).__init__()

        self._input_path = input_path
        self._proc_id = proc_id
        
        self._outqu = outqu
        self._pause_qu = pause_qu

        self._Task = Task
        self._task_n = task_n

        self._constants = constants
        self._task_size = constants.task_size
        self._reader_n = constants.reader_n
        self._temp_dir = constants.temp_dir
        self._debug = constants.debug

    #Find data pertaining to assocd and all reads 
    #and divide pertaining to the chromosome that it is aligned to
    def run(self):
        parent_bam_mem_obj = ParentAlignmentFile(self._input_path)
        parent_bam = pysam.AlignmentFile(self._input_path,"rb")
        parent_iter = self.__get_parent_iter__(parent_bam)
        parent_generator = self.__bam_generator__(parent_iter)
        wait_for_pause = self.__wait_for_pause__
        if self._debug:
            parent_generator = self.__debug_generator__(parent_iter)

        task_qu = Queue()

        tasks = [self._Task(parent_bam=parent_bam_mem_obj,
                      inqu=task_qu,
                      outqu=self._outqu,
                      task_size=self._task_size,
                      constants=self._constants) \
                        for i in xrange(self._task_n)]

        for task in tasks:
            task.start()

        for i,command in enumerate(parent_generator):
            #This pause adds stability to parabam
            #I think this is because it spaces calls to disk
            time.sleep(.0001) 
            wait_for_pause()
            task_qu.put(parent_bam.tell())

        for n in xrange(self._task_n+1):
            task_qu.put(DestroyPackage())

        time.sleep(3)
        parent_bam.close()
        task_qu.close()
        return

    def __get_parent_iter__(self,parent_bam):
        if not self._constants.fetch_region:
            return parent_bam.fetch(until_eof=True)
        else:
            return parent_bam.fetch(region=self._constants.fetch_region)

    def __bam_generator__(self,parent_iter):
        cdef int proc_id = self._proc_id
        cdef int reader_n = self._reader_n
        cdef int iterations = 0
        cdef int task_size = self._task_size

        while True:            
            try:
                if iterations % reader_n == proc_id:
                    yield True
                for x in xrange(task_size):
                    parent_iter.next()
                iterations += 1
            except StopIteration:
                break
        return

    def __debug_generator__(self,parent_iter):
        cdef int proc_id = self._proc_id
        cdef int reader_n = self._reader_n
        cdef int iterations = 0
        cdef int task_size = self._task_size

        while True:            
            try:
                if iterations % reader_n == proc_id:
                    yield True
                for x in xrange(task_size):
                    parent_iter.next()
                iterations += 1

                if iterations == 25:
                    break

            except StopIteration:
                break
        return

    def __query_pause_qu__(self,bypass,pause_qu):
        #This fairly obfuscated code is to fix
        #a race condition when a pause can
        #be missed due to waiting on a previous pause
        try:
            if bypass:
                pause = True
            else:
                pause = pause_qu.get(False)

            if pause:
                while True:
                    try:
                        pause = pause_qu.get(False)
                        if not pause:
                            break
                    except Queue2.Empty:
                        time.sleep(1)
        except Queue2.Empty:
            pass

        last = False  
        while True:
            try:
                pause = pause_qu.get(False)
                last = pause
            except Queue2.Empty:
                break
        return last

    def __wait_for_pause__(self):
        pause_qu = self._pause_qu
        pause = self.__query_pause_qu__(False,pause_qu)
        if pause:
            while True:
                pause = self.__query_pause_qu__(True,pause_qu)
                if not pause:
                    break
        
class Leviathon(object):
    #Leviathon takes objects of file_readers and handlers and
    #chains them together.
    def __init__(self,object constants,dict handler_bundle,
                 list handler_order,list queue_names,
                 int update,object Task):

        self._constants = constants
        self._handler_bundle  = handler_bundle
        self._handler_order = handler_order
        self._update = update
        self._queue_names = queue_names
        self._Task = Task
    
    def run(self,input_path,output_paths):
        parent = ParentAlignmentFile(input_path)

        default_qus = self.__create_queues__(self._queue_names) 
        pause_qus = self.__create_pause_qus__(self._constants.reader_n)


        handlers_objects,handler_inqus = self.__create_handlers__(self._handler_order,
                                                          self._handler_bundle,self._constants,
                                                          default_qus,parent,output_paths,pause_qus)
        handlers = self.__get_handlers__(handlers_objects)

        task_n = self.__get_task_n__(self._constants,handlers)
        file_reader_bundles = self.__get_file_reader_bundles__(default_qus,parent,
                                                               task_n,self._constants,pause_qus)
        file_readers = self.__get_file_readers__(file_reader_bundles)

        #Start file_readers
        for file_reader in file_readers:
            file_reader.start()
            #time.sleep(2)

        #Start handlers:
        for handler in handlers:
            handler.start()

        #Wait for file_readers to finish
        for file_reader in file_readers:
            file_reader.join()

        #Destory handlers
        for handler,queue in izip(handlers,handler_inqus):
            queue.put(DestroyPackage())
            handler.join()

        del default_qus
        del file_reader_bundles
        del file_readers
        del pause_qus
        del handler_inqus
        del handlers_objects
        del handlers

        gc.collect()


    def __create_pause_qus__(self,reader_n):
        return [Queue() for i in xrange(reader_n)]

    def __get_task_n__(self,constants,handlers):
        task_n = (constants.total_procs - len(handlers) - constants.reader_n) / constants.reader_n
        if task_n > 0:
            return task_n
        else:
            return 1

    def __get_file_readers__(self,file_reader_bundles):
        file_readers = []
        for bundle in file_reader_bundles: 
            file_reader = FileReader(**bundle)
            file_readers.append(file_reader)
        return file_readers

    def __proc_id_generator__(self,reader_n):
        for i in reversed(xrange(reader_n)):
            yield i

    def __get_file_reader_bundles__(self,default_qus,parent,task_n,constants,pause_qus):
        bundles = []

        for proc_id,pause in izip(self.__proc_id_generator__(constants.reader_n),pause_qus):
            current_bundle = {"input_path" :parent.filename,
                              "proc_id" :proc_id,
                              "task_n" :task_n,
                              "outqu" :default_qus["main"],
                              "constants" :constants,
                              "Task" :self._Task,
                              "pause_qu" :pause}

            bundles.append(current_bundle)

        return bundles
 
    def __create_queues__(self,queue_names):
        queues = {}
        for name in queue_names:
            queues[name] = Queue()
        return queues

    def __create_handlers__(self,handler_order,handler_bundle,constants,
                            queues,parent_bam,output_paths,pause_qus):
        handlers = []
        handler_inqus = []

        for handler_class in handler_order:
            handler_args = dict(handler_bundle[handler_class])

            handler_args["parent_bam"] = parent_bam
            handler_args["output_paths"] = output_paths
            handler_args["constants"] = constants

            #replace placeholder with queues
            handler_args["pause_qus"] = pause_qus
            handler_args["inqu"] = queues[handler_args["inqu"]]
            handler_args["out_qu_dict"] = dict(\
                    [(name,queues[name]) for name in handler_args["out_qu_dict"] ])

            handler_inqus.append(handler_args["inqu"])
            handlers.append(handler_class(**handler_args))

        return handlers,handler_inqus

    def __get_handlers__(self,handlers):
        handler_processes = []
        for handler in handlers:
            hpr = Process(target=handler.listen,args=(self._update,))
            handler_processes.append(hpr)
        return handler_processes

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
            ,help="The maximum amount of processes you wish parabam to use. This should"\
                  "be less than or equal to the amount of processor cores in your machine.")
        parser.add_argument('-s',type=int,nargs='?',default=250000
            ,help="The amount of reads considered by each distributed task.")
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

class Constants(object):
    
    def __init__(self,temp_dir,verbose,task_size,total_procs,fetch_region,**kwargs):
        #TODO: Update with real required constants values
        self.temp_dir = temp_dir
        self.verbose = verbose
        self.task_size = task_size
        self.total_procs = total_procs
        self.fetch_region = fetch_region

        for key, val in kwargs.items():
            setattr(self,key,val)

    def add(self,key,val):
        setattr(self,key,val)

class Package(object):
    def __init__(self,results):
        self.results = results

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

        if has_index:
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
