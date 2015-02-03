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

    def __init__(self,object task_set,object outqu,int curproc,object destroy,object const,str parent_class):
        multiprocessing.Process.__init__(self)
        self._task_set = task_set
        self._outqu = outqu
        self._curproc = curproc
        self._temp_dir = const.temp_dir
        self._destroy = destroy
        self._parent_class = parent_class
        self.const = const

    def run(self):
        cdef dict results = self.__generate_results__()
        results["total"] = len(self._task_set)

        self._outqu.put(CorePackage(name=self.name,
                                results=results,
                                destroy=self._destroy,
                                curproc=self._curproc,
                                parent_class=self._parent_class))
        #Trying to control mem useage
        del self._task_set
        gc.collect()

    def __get_temp_path__(self,typ,ext=".bam"):
        #self.pid ensures that the temp names are unique.
        return "%s/%s_%s_parabam_temp%s" % (self._temp_dir,typ,self.pid,ext)

    def __get_mate__(self,master,alig):
        try:
            start = master.tell()
            mate = master.mate(alig)
            master.seek(start)
        except ValueError:
            mate = None
        return mate

    #Generate a dictionary of results using self._task_set
    @abstractmethod
    def __generate_results__(self):
        #Should always return a dictionary
        pass

cdef class Handler:

    def __init__(self,object inqu,object const,object report=True,int destroy_limit=1):
        self._inqu = inqu
        self._destroy_limit = destroy_limit
        self._report = report
        self.const = const

        self._verbose = const.verbose
        self._output_paths = const.output_paths

        self._periodic_interval = 10
        self._destroy_count = 0

        self._stats = {}

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
        destroy = False

        cdef int iterations = 0
        cdef int start_time = time.time()
        cdef int dealt  = 0
        cdef int curproc = 0

        update_output = self._update_output #speedup alias
        periodic_interval = self._periodic_interval #speedup alias

        while not destroy:
            iterations += 1
            #Listen for a process coming in...
            try:
                new_package = self._inqu.get(False)
                if new_package.destroy:
                    #destroy limit allows for multiple processors
                    self._destroy_count += 1
                    if self._destroy_count == self._destroy_limit:
                        destroy = True

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

    def __init__(self,object outqu,object const,object TaskClass,list task_args,object debug=False):

        self._debug = debug

        self._chunk = const.chunk
        self._proc = const.proc

        self._master_file_path = const.master_file_path
        self._verbose = const.verbose
        self._temp_dir = const.temp_dir

        self._outqu = outqu
        self.const = const

        #Class which we wish to run on each read
        self._TaskClass = TaskClass
        self._task_args = task_args 

        self._active_tasks = []
        self._active_count = 0
        self._prev = 0

    #Find data pertaining to assocd and all reads 
    #and divide pertaining to the chromosome that it is aligned to
    def run(self,int update):
    
        self.__pre_processor__(self._master_file_path)
        
        cdef object master_bam = self.__get_master_bam__(self._master_file_path)
        cdef list collection = []
        cdef int collection_count = 0

        #function alias, for optimisation
        add_to_collection = self.__add_to_collection__
        start_task = self.__start_task__
        wait_for_tasks = self.__wait_for_tasks__

        #Insert a debug iterator into the processor, this workaround doesn't
        #slow processing when not in debug mode.
        bam_iterator = self.__get_next_alig__(master_bam) if not self._debug \
                            else self.__get_next_alig_debug__(master_bam)

        for i,alig in enumerate(bam_iterator):

            add_to_collection(master_bam,alig,collection)
            collection_count += 1

            if collection_count == self._chunk:
                wait_for_tasks(self._active_tasks,self._proc) # -2 for proc and handler 
                start_task(collection)      
                del collection
                collection = [] # already running
                collection_count = 0

        wait_for_tasks(self._active_tasks,0)
        start_task(collection,destroy=True)
        self.__end_processing__(master_bam)

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
        args = [collection,self._outqu,self._active_count,destroy,self.const,self.__class__.__name__]
        args.extend(self._task_args)
        task = self._TaskClass(*args)
        task.start()
        self._active_tasks.append(task)

    def __get_master_bam__(self,master_file_path):
        return pysam.Samfile(master_file_path,"rb") #Open telbam for analysis

    #Should be over written if we don't want data from BAM file.
    def __get_next_alig__(self,master_bam):
        for alig in master_bam.fetch(until_eof=True):
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
        #del self._active_tasks
        #gc.collect()

    def __pre_processor__(self,master_bam):
        pass

    def __add_to_collection__(self,master,item,collection):
        pass

class Leviathon(object):
    #Leviathon takes objects of processors and handlers and
    #chains them together.
    def __init__(self,list processors,list handlers,int update=10000000):
        self._processors = processors 
        self._handlers  = handlers
        self._update = update
    
    def run(self):
        procs = []

        for p in self._processors:
            ppr = multiprocessing.Process(target=p.run,args=(self._update,))
            ppr.start()
            procs.append(ppr)

        for h in self._handlers:
            hpr = multiprocessing.Process(target=h.listen,args=(self._update,))
            hpr.start()
            procs.append(hpr)

        procs[-1].join() #Wait on the last handler that we started.

        for proc in procs:
            proc.join()
            proc.terminate()

        del self._processors
        del self._handlers
        del procs

#Provides a conveinant way for providing an Interface to parabam
#programs. Includes default command_args and framework for 
#command-line and programatic invocation. 
class Interface(object):

    __metaclass__ = ABCMeta

    def __init__(self,temp_dir,exe_dir):
        self._temp_dir = temp_dir
        self._exe_dir = exe_dir

    def __introduce__(self,name):
        intro =  "%s has started. Start Time: " % (name,)\
            + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        underline = "-" * len(intro)
        print intro
        print underline

    def __goodbye__(self,name):
        print "%s has finished. End Time: " % (name,)\
            + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')

    def __get_group__(self,bams,names,multi):
        if multi == 0:
            multi = 1
        for i in xrange(0,len(bams),multi):
            yield (bams[i:i+multi],names[i:i+multi])

    def __get_basename__(self,path):
        base = os.path.basename(path)
        if "." in base:
            base = base.rpartition(".")[0]
        return base

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

    def default_parser(self):

        parser = argparse.ArgumentParser(conflict_handler='resolve',
                    formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-p',type=int,nargs='?',default=4
            ,help="The maximum amount of tasks run concurrently. This number should\n"\
            "be less than or equal to the amount of cores on your machine [Default: 4]")
        parser.add_argument('-c',type=int,nargs='?',default=100000
            ,help="How many reads each parallel task is given. A higher number\n"\
            "will mean faster analysis but also a greater burden on memory [Default:100000]")
        parser.add_argument('-v',action="store_true",default=False
            ,help='If set the program will output more information to terminal')
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

class UserInterface(Interface):

    __metaclass__ = ABCMeta

    def __init__(self,temp_dir,exe_dir):
        super(UserInterface,self).__init__(temp_dir,exe_dir)

    def __get_module_and_vitals__(self,code_path):
        if os.getcwd not in sys.path:
            sys.path.append(os.getcwd())
        try:
            module = __import__(code_path, fromlist=[''])

        except ImportError:
            sys.stderr.write("[ERROR] parabam can't find user module. Ensure code is in current working directory\n")
            raise SystemExit

        user_engine = module.engine
        user_constants = {}
        if hasattr(module,"set_constants"):
            module.set_constants(user_constants)

        return module,user_engine,user_constants

    def default_parser(self):
        parser = super(UserInterface,self).default_parser()

        parser.add_argument('--instruc','-i',metavar='INSTRUCTION',required=True
            ,help='The instruction file, written in python, that we wish'\
            'to carry out on the input BAM.')
        parser.add_argument('--input','-b',metavar='INPUT', nargs='+',required=True
            ,help='The file(s) we wish to operate on. Multipe entries should be separated by a single space')
        return parser
        
    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def get_parser(self):
        pass

class Const(object):
    
    def __init__(self,output_paths,temp_dir,master_file_path,verbose,chunk,proc,**kwargs):
        self.output_paths = output_paths
        self.temp_dir = temp_dir
        self.master_file_path = master_file_path
        self.verbose = verbose
        self.chunk = chunk
        self.proc = proc

        for key, val in kwargs.items():
            setattr(self,key,val)

    def add(self,key,val):
        setattr(self,key,val)

class Package(object):
    def __init__(self,name,results,destroy):
        self.name = name
        self.results = results
        self.destroy = destroy

class CorePackage(Package):
    def __init__(self,name,results,destroy,curproc,parent_class):
        super(CorePackage,self).__init__(name,results,destroy)
        self.curproc = curproc
        self.parent_class = parent_class

#And they all lived happily ever after...
