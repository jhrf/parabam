import os
import time
import datetime
import sys
import Queue as Queue2
import gc
import shutil
import argparse
import tempfile
import signal
import traceback

import pysam

from itertools import izip
from multiprocessing import Process,Queue,freeze_support
from abc import ABCMeta, abstractmethod

class CmdLineInterface(object):
    def __init__(self,program_name = "parabam"):
        self.program_name = program_name

    def keyboard_handler(self,sig, frame):#Catch keyboard interupt and end processors
        sys.stdout.write("\r[ERROR] %s interrupted by user" \
                                         % (self.program_name))
        sys.exit(0)

    def die_gracefully(self, interface):
        if interface is not None:
            temp_dir = interface.get_temp_dir_path()
            if not temp_dir == "." or not temp_dir == "" or not temp_dir == "./": 
                interface.interface_exit()

    def handle(self,command_map,help_text):
        freeze_support()

        signal.signal(signal.SIGINT, self.keyboard_handler)     
        print_help_text = False

        if len(sys.argv) >= 2 and sys.argv[1] in command_map:
            
            command = sys.argv[1]
            #Remove command from arguments.
            sys.argv = sys.argv[1:] 
            interface = None
            try:
                #Load the command using the command line
                interface = command_map[command](cmd_run=True)
                interface.run_cmd()
                
            except SystemExit:
                print " "
                print "[Status] %s is quitting gracefully\n" \
                                                % (self.program_name,)
            except BaseException as exception:
                print " "
                print "[Error] %s stopped unexpecedtly, sorry!" \
                                                % (self.program_name,)

                traceback.print_exception(*sys.exc_info())
            finally:
                self.die_gracefully(interface)
        
        else: 
            if len(sys.argv) >= 2 and sys.argv[1] not in command_map:
                print "\nUnrecognised command: `%s`." % (sys.argv[1])
                print "Refer to manual below for a list of valid commands\n"
            print help_text

cdef class Handler:

    def __init__(self,
                  object parent_bam, 
                  object output_paths,
                  object inqu,
                  object constants,
                  object pause_qus,
                  dict out_qu_dict,
                  object report=True):

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

        self._processing = True
        self._destroy = False
        self._finished = False
        self._last_update = False

        if constants.verbose == 1:
            self._verbose = True
            self._update_output = self.__level_1_output__
        elif constants.verbose == 2:
            #In the case verbose is simply "True" or "level 2"
            self._verbose = True
            self._update_output = self.__level_2_output__
        else:#catching False and -v0
            self._verbose = False
            self._update_output = self.__level_1_output__

    def __standard_output__(self,out_str):
        sys.stdout.write(out_str + "\n")
        sys.stdout.flush()

    def __level_1_output__(self,out_str):
        total_procd = self.__total_reads__()
        time = out_str.partition("Time: ")[2]
        sys.stdout.write("\t- Reads processed: %d Time: %s\n" \
                                                % (total_procd,time))
        sys.stdout.flush()

    def __level_2_output__(self,outstr):
        sys.stdout.write("\r" + outstr)
        sys.stdout.flush()

    #Must be overwritten if stats architecture is modififed
    def __total_reads__(self):
        if self._stats == {}:
            return 0
        return self._stats["Total"]

    def listen(self,update_interval):
        destroy = self._destroy

        cdef int iterations = 0
        cdef int start_time = time.time()
        cdef int dealt  = 0

        periodic_interval = self._periodic_interval #speedup alias
        finished = self.__is_finished__

        output_logic = self.__silence__
        if self._constants.verbose and self._report:
            output_logic = self.__verbose_output_logic__

        while not finished():
            iterations += 1
            #Listen for a process coming in...
            try:
                packages = []
                for i in xrange(20):
                    packages.append(self._inqu.get(False))
            except Queue2.Empty:
                #Queue empty. Continue with loop
                #time.sleep(.5)
                time.sleep(.01)

            for package in packages:
                if not package.results == {}:#If results are present...
                    self.__new_package_action__(package) #Handle the results
                    dealt += 1
                if type(package) == DestroyPackage:
                    self._destroy = True
                elif type(package) == EndProcPackage:
                    self._processing = False
            
            del packages
            if iterations % periodic_interval == 0: 
                self.__periodic_action__(iterations)

            output_logic(start_time,iterations,update_interval)

        self._inqu.close()
        self.__handler_exit__()

    def __verbose_output_logic__(self,start_time,iterations,update_interval):

        if iterations % update_interval == 0 and self._processing:
            outstr = self.__format_update__(start_time,
                                            iterations,
                                            update_interval)
            self._update_output(outstr)
        elif not self._processing and not self._last_update:
            if self._constants.verbose == 2:
                outstr = self.__format_update__(start_time,
                                                iterations,
                                                update_interval)
                self._update_output(outstr)
                sys.stdout.write("\n")
                sys.stdout.flush()

            self._last_update = True

    def __silence__(self,*args):
        pass

    def __is_finished__(self):
        if not self._destroy or not self._finished:
            return False
        return True

    def __format_update__(self,start_time,iterations,update_interval):
        first_char = self.__get_first_char__(iterations,update_interval)
        stats = "\t%s " % (first_char,)
        for name,value in self._stats.items():
            stats += "%s:%d " % (name,value)
        return "\r%s| Time: %d " %\
                (stats,self.__secs_since__(start_time),)

    def __get_first_char__(self,iterations,update_interval):
        first_chars = ["-","+","*","+"]
        return first_chars[(iterations/update_interval) % len(first_chars)]

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

    def __init__(self, parent_bam, inqu, outqu, statusqu, task_size, constants):
        super(Task,self).__init__()
        self._parent_bam = parent_bam
        self._parent_path = self._parent_bam.filename
        self._header = self._parent_bam.header

        self._outqu = outqu
        self._inqu = inqu
        self._statusqu = statusqu

        self._task_size = task_size
        self._temp_dir = constants.temp_dir
        self._constants = constants 
        
        self._dealt = 0

    def run(self):

        bamfile = pysam.AlignmentFile(self._parent_path,"rb")
        iterator = bamfile.fetch(until_eof=True)
        generate_results = self.__generate_results__

        while True:
            try:
                package,sequence_id = self._inqu.get(False)
                self._idle_sent = False
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

                self._outqu.put(Package(results=results,sequence_id=sequence_id))
                #tell filereader a batch of five jobs
                self._statusqu.put(1) 

            except Queue2.Empty:
                time.sleep(5)
            except StopIteration:
                results = self.__handle_stop_iteration__(seek,bamfile,iterator)
                self._dealt += 1 #Copy paste code here.
                time.sleep(0.005) #And here
                self.__post_run_routine__()
                self._outqu.put(Package(results=results,sequence_id=sequence_id))
        return

    def __generate_results__(self,iterator,**kwargs):
        self.__pre_run_routine__(iterator)
        self.__process_task_set__(iterator)
        results = self.__get_results__()
        results["Total"] = self._task_size
        self.__post_run_routine__()
        return results

    def __handle_stop_iteration__(self,seek,bamfile,iterator):
        results = self.__get_results__()
         
        cdef int count = 0 
        bamfile.seek(seek)
        try: #Count the amount of reads until end of file
            for read in iterator:
                count += 1
        except StopIteration:
            pass

        results["Total"] = count
        return results

    def __get_temp_object__(self,path):
        return pysam.AlignmentFile(path,"wb",header=self._parent_bam.header)

    def __get_temp_path__(self,identity):
        file_name = "%s_%d_%d_%s" %\
            (identity,self.pid,self._dealt, 
             os.path.split(self._parent_bam.filename)[1])

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
                 object pause_qu,object inqu = None):
        
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

        if not inqu:
            self._inqu = Queue()
        else:
            self._inqu = inqu
        self._active_jobs = 0
        self._active_jobs_thresh = ((self._task_n * 50) *\
                                        ((1000/self._task_size)+1))

    #Find data pertaining to assocd and all reads 
    #and divide pertaining to the chromosome that it is aligned to
    def run(self):

        parent_bam_mem_obj = ParentAlignmentFile(self._input_path)
        parent_bam = pysam.AlignmentFile(self._input_path,"rb")
        parent_iter = self.__get_parent_iter__(parent_bam)
        parent_generator = self.__bam_generator__(parent_iter)
        wait_for_pause = self.__wait_for_pause__
        check_inqu = self.__check_inqu__


        if self._debug:
            parent_generator = self.__debug_generator__(parent_iter)

        task_qu = Queue()

        tasks = [self._Task(parent_bam=parent_bam_mem_obj,
                      inqu=task_qu,
                      outqu=self._outqu,
                      statusqu = self._inqu,
                      task_size=self._task_size,
                      constants=self._constants) \
                        for i in xrange(self._task_n)]

        for task in tasks:
            task.start()

        for i,command in enumerate(parent_generator):
            wait_for_pause()
            task_qu.put( (parent_bam.tell(),command) )
            self._active_jobs += 1
            if i % 10 == 0:
                check_inqu()

        for n in xrange(self._task_n+1):
            task_qu.put( (DestroyPackage(),-1) )

        time.sleep(2)
        parent_bam.close()
        task_qu.close()
        return

    def __check_inqu__(self):

        while True:
            try:
                self._active_jobs -= self._inqu.get(False)
            except Queue2.Empty:
                break

        if self._active_jobs > self._active_jobs_thresh:
            while True:
                try: #Block until active jobs is less than limit
                    self._active_jobs -= self._inqu.get(False)
                    if self._active_jobs <= self._task_n*2:
                        break
                except Queue2.Empty:
                    self.__wait_for_pause__()
                    time.sleep(1)

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
                    yield iterations
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
                    yield iterations
                for x in xrange(task_size):
                    parent_iter.next()
                iterations += 1

                if iterations == 25:
                    break

            except StopIteration:
                break
        return

    def __send_ack__(self,qu):
        qu.put(2)
        time.sleep(1)

    def __wait_for_pause__(self):
        try: #pre clause for speed
            pause = self._pause_qu.get(False)
        except Queue2.Empty:
            return

        pause_qu = self._pause_qu
        while True:
            try:#get most recent signal
                attempt = pause_qu.get(False)
                pause = attempt
            except Queue2.Empty:
                break
        
        if pause == 1:
            self.__send_ack__(pause_qu)
            while True:
                try:
                    pause = pause_qu.get(False) 
                    if pause == 0: #Unpause signal recieved
                        self.__send_ack__(pause_qu)
                        return
                except Queue2.Empty:
                    time.sleep(.5)
 
class Leviathan(object):
    #Leviathon takes objects of file_readers and handlers and
    #chains them together.
    def __init__(self,object constants,dict handler_bundle,
                 list sequence_id,list queue_names,
                 int update,object Task,object FileReaderClass):

        self._constants = constants
        self._handler_bundle  = handler_bundle
        self.sequence_id = sequence_id
        self._update = update
        self._queue_names = queue_names
        self._Task = Task
        self._FileReaderClass = FileReaderClass
    
    def run(self,input_path,output_paths):
        parent = ParentAlignmentFile(input_path)

        default_qus = self.__create_queues__(self._queue_names) 
        pause_qus = self.__create_pause_qus__(self._constants.reader_n)

        handlers_objects,handler_inqus = \
                        self.__create_handlers__(self.sequence_id,
                                                self._handler_bundle,
                                                self._constants,
                                                default_qus,
                                                parent,
                                                output_paths,
                                                pause_qus)                        
        handlers = self.__get_handlers__(handlers_objects)

        task_n_list = self.__get_task_n__(self._constants,handlers)
        file_reader_bundles = self.__get_file_reader_bundles__(default_qus,
                                                                parent,
                                                                task_n_list,
                                                                self._constants,
                                                                pause_qus)
        file_readers = self.__get_file_readers__(file_reader_bundles)

        #Start file_readers
        for file_reader in file_readers:
            file_reader.start()

        #Start handlers:
        for handler in handlers:
            handler.start()

        #Wait for file_readers to finish
        for file_reader in file_readers:
            file_reader.join()

        #Inform handlers that processing has finished
        for handler,queue in izip(handlers,handler_inqus):
            queue.put(EndProcPackage())

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
        total_procs = constants.total_procs
        reader_n = constants.reader_n
        task_n_list = [total_procs / reader_n] * reader_n

        if not total_procs % reader_n == 0:
            for i,overflow in zip(xrange(len(task_n_list)),
                                  [1] * (total_procs % reader_n)):
                task_n_list[i] += overflow

        if any([ n == 0 for n in task_n_list]):
            task_n_list = [1] * reader_n

        return task_n_list 

    def __get_file_readers__(self,file_reader_bundles):
        file_readers = []
        for bundle in file_reader_bundles: 
            file_reader = self._FileReaderClass(**bundle)
            file_readers.append(file_reader)
        return file_readers

    def __proc_id_generator__(self,reader_n):
        for i in reversed(xrange(reader_n)):
            yield i

    def __get_file_reader_bundles__(self, 
                                     default_qus, 
                                     parent, 
                                     task_n_list, 
                                     constants, 
                                     pause_qus):
        bundles = []

        for proc_id,pause,task_n in izip(
                                self.__proc_id_generator__(constants.reader_n),
                                pause_qus,
                                task_n_list):

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

    def __create_handlers__(self,sequence_id,handler_bundle,constants,
                            queues,parent_bam,output_paths,pause_qus):
        handlers = []
        handler_inqus = []

        for handler_class in sequence_id:
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

class Interface(object):
    '''Provides a conveinant way for providing an Interface to parabam
    programs. Includes default command_args and framework for
    command-line and programatic invocation.'''

    __metaclass__ = ABCMeta

    def __init__(self,
                  instance_name = "parabam",
                  temp_dir = None,
                  keep_in_temp = False,
                  task_size = 250000,
                  total_procs = 8,
                  reader_n = 2,
                  verbose = False,
                  announce = False,
                  cmd_run = False):

        self.instance_name = instance_name
        self.header_line = "-" * 77

        self.cmd_run = cmd_run
        if cmd_run:
            self.__cmd_args_to_class_vars__()

        else:
            self.task_size = task_size
            self.total_procs = total_procs
            self.reader_n = reader_n
            self.verbose = verbose
            self.keep_in_temp = keep_in_temp
            if not verbose:
                self.announce = False
            else:
                self.announce = announce

        self.temp_dir = temp_dir
        self.__temp_dir_instalise__()

    def get_temp_dir_path(self):
        return self.temp_dir

    def __temp_dir_instalise__(self):
        if self.temp_dir is None:
            self.temp_dir = self.__get_unique_tempdir__()
            self.control_temp_dir = True
        else:
            self.temp_dir = self.temp_dir
            self.control_temp_dir = False

    def __get_unique_tempdir__(self):
        sanitised_name = self.instance_name.replace(" ","_")
        prefix = "%s-" % (sanitised_name,)
        return tempfile.mkdtemp(prefix=prefix,dir=".")

    def __cmd_args_to_class_vars__(self):
        parser = self.get_parser()
        cmd_args = parser.parse_args()

        self.task_size = cmd_args.s
        self.total_procs = cmd_args.p
        self.verbose = cmd_args.v
        self.reader_n = cmd_args.f
    
        self.keep_in_temp = False
        self.announce = int(self.verbose) > 0

        self.cmd_args = cmd_args

    def __introduce__(self):
        if self.verbose and self.announce:
            intro =  "%s has started. Start Time: %s" \
                        % (self.instance_name,self.__get_date_time__())
            underline = ("-" * len(intro))

            sys.stdout.write(intro+"\n")
            sys.stdout.write(underline+"\n\n")
            sys.stdout.flush()
    
    def __get_date_time__(self):
        return datetime.datetime.fromtimestamp(
                            time.time()).strftime('%Y-%m-%d %H:%M:%S')

    def __goodbye__(self):
        if self.verbose and self.announce:
            outro = "%s has finished.  End Time: %s" \
                        % (self.instance_name,self.__get_date_time__())

            sys.stdout.write(("-" * len(outro)) + "\n")
            sys.stdout.write(outro + "\n")
            sys.stdout.flush()

    def interface_exit(self):
        if self.control_temp_dir and \
                not self.keep_in_temp and\
                not self.temp_dir is None:

            for path in os.listdir(self.temp_dir):
                remove_path = os.path.join(self.temp_dir,path)
                os.remove(remove_path)
            shutil.rmtree(self.temp_dir)
            self.temp_dir = None

    def default_parser(self):

        parser = ParabamParser(conflict_handler='resolve',
                    formatter_class=argparse.RawTextHelpFormatter)

        parser.add_argument('-p',type=int,nargs='?',default=4
            ,help=('The maximum amount of processes you wish\n' 
                   '%s to use. This should be less \n'
                   'than or equal to the amount of processor\n'
                   'cores in your machine [Default: 4].') % self.instance_name)
        parser.add_argument('-s',type=int,nargs='?',default=250000
            ,help="The amount of reads considered by each\n"\
                    "distributed task. [Default: 250000]")
        parser.add_argument('-f',type=int,metavar="READERS",nargs='?',default=1,
            help="The amount of open connections to the file being read.\n"\
                    "Conventional hard drives perform best with \n"\
                    "the default of 1. [Default: 1]")
        parser.add_argument('-v', choices=[0,1,2],default=0,type=int,
            help=("The amount of information output by the program:\n"
                 "\t0: No output [Default]\n"
                 "\t1: Total Reads Processed\n"
                 "\t2: Detailed output"))
        
        return parser

    @abstractmethod
    def run_cmd(self):
        #This is usualy just a function that
        #takes an argparse parser and turns 
        #passes the functions to the run function
        pass

    @abstractmethod
    def run(self):
        #The main run method. Should make a call to `interface_exit`
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
    
    def __init__(self, **kwargs):

        for key, val in kwargs.items():
            setattr(self,key,val)

    def add(self,key,val):
        setattr(self,key,val)

class Package(object):
    def __init__(self,results,sequence_id=-1):
        self.results = results
        self.sequence_id = sequence_id

class DestroyPackage(Package):
    def __init__(self):
        super(DestroyPackage,self).__init__(results={})
        self.destroy = True

class EndProcPackage(Package):
    def __init__(self):
        super(EndProcPackage,self).__init__(results={})
        self.destroy = False

class ParentAlignmentFile(object):
    
    def __init__(self,path):
        has_index = os.path.exists(os.path.join("%s%s" % (path,".bai")))

        parent = pysam.AlignmentFile(path,"rb")
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
