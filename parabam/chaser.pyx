#Once upon a time
import pysam
import time
import sys
import os
import gc
import parabam
import random

import Queue as Queue2
import numpy as np

from parabam.core import DestroyPackage
from itertools import izip
from collections import Counter,namedtuple
from multiprocessing import Queue,Process

class ChaserResults(parabam.core.Package):
    def __init__(self,results,chaser_type):
        super(ChaserResults,self).__init__(results)
        self.chaser_type=chaser_type

class PrimaryResults(ChaserResults):
    def __init__(self,results,chaser_type):
        super(PrimaryResults,self).__init__(results,chaser_type)

class MatchMakerResults(ChaserResults):
    def __init__(self,loner_type,results,chaser_type,rescued,level):
        super(MatchMakerResults,self).__init__(results,chaser_type)
        self.rescued = rescued
        self.loner_type = loner_type
        self.level = level

class JobPackage(object):
    def __init__(self,paths,loner_type,level):
        self.paths = paths
        self.loner_type = loner_type 
        self.level = level

class Handler(parabam.core.Handler):

    def __init__(self,object parent_bam, object output_paths,
                      object inqu,object constants,object pause_qus,
                      dict out_qu_dict,object Task):

        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                             inqu=inqu,constants=constants,pause_qus=pause_qus,
                             out_qu_dict=out_qu_dict,report = False)
        
        class EngineTask(Task):
            def __init__(self,parent_bam,constants):
                super(EngineTask,self).__init__(parent_bam=parent_bam,
                                                    inqu=None,
                                                    outqu=None,
                                                    statusqu=None,
                                                    task_size=0,
                                                    constants=constants)

            def handle_pair_dict(self,pairs,unique):
                self.unique = unique
                results = self.__generate_results__(pairs)
                results["total"] = 0
                self._dealt += 1
                time.sleep(.01)
                return results

            def __process_task_set__(self,iterator):
                engine = self._engine
                parent_bam = self._parent_bam
                handle_output = self.__handle_engine_output__
                user_constants = self._user_constants

                for read_name,pair in iterator.items():
                    engine_output = engine(pair,user_constants,parent_bam)
                    handle_output(engine_output,pair)

            def __get_temp_path__(self,identity):
                file_name = "%s_%d_%d_%s" %\
                    (identity,self.unique,self._dealt, os.path.split(self._parent_bam.filename)[1])
                return os.path.join(self._temp_dir,file_name)

        self._loner_pyramid = self.__instalise_loner_pyramid__()
        self._pyramid_idle_counts = Counter()
        self._loner_purgatory = {}

        self._periodic_interval = 5
        self._total_loners = 0
        self._rescued = Counter()
        self._prev_rescued = 0
        self._stale_count = 0

        self._child_pack = {"parent_bam":parent_bam,
                            "queue":out_qu_dict["main"],
                            "chrom_bins":self.__get_chrom_bins__(parent_bam) ,
                            "chaser_task":EngineTask(parent_bam,self._constants)}

        self._chaser_qu = Queue()
        self._chaser_tasks = self.__create_chaser_tasks__(self._constants.total_procs)
        self._pending_jobs = 0
        self._max_jobs = (self._constants.total_procs * 3)
        self._file_readers_paused = False

        self._post_destroy_count = 0

        self._primary_store = self.__instalise_primary_store__()

    def __create_chaser_tasks__(self,total_tasks):
        tasks = []
        for i in xrange(total_tasks):
            task = ChaserTask(constants=self._constants,
                              inqu=self._chaser_qu,
                              outqu=self._inqu,
                              max_reads=7000000,
                              child_pack=self._child_pack)
            task.start()
            tasks.append(task)
        return tasks

    def __get_chrom_bins__(self,parent_bam):
        references = parent_bam.references
        lengths = parent_bam.lengths

        chrom_bins = {}
        chrom_info = zip(xrange(len(references)),references,lengths)
        chrom_info = sorted(chrom_info, key=lambda tup: tup[2])
        group = []
        size = 0
        uid = 0
        for i,name,length in chrom_info:
           group.append( (i,name,length,) )
           size += length

           if size > 150000000:
                for i,name,length in group:
                    chrom_bins[i] = "b%d" % (uid,)
                group = []
                size = 0
                uid += 1
        return chrom_bins
           
    def __instalise_primary_store__(self):
        return []

    def __instalise_loner_pyramid__(self):
        return {}

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

        return len(active_tasks)

    def __new_package_action__(self,new_package,**kwargs):        
        chaser_type = new_package.chaser_type
        if chaser_type == "origin":
            self.__handle_origin_task__(new_package)
        elif chaser_type == "match_maker":
            self.__handle_match_maker__(new_package)
            self._pending_jobs -= 1
        elif chaser_type == "primary":
            self.__handle_primary_task__(new_package)
            self._pending_jobs -= 1
        else:
            print "Error!"

    def __handle_origin_task__(self,new_package):
        for loner_count,path in new_package.results:
            if loner_count == 0:
                os.remove(path)
            else:
                self._total_loners += loner_count
                self._primary_store.append(path) 

    def __handle_match_maker__(self,new_package):
        loner_type = new_package.loner_type
        level = new_package.level

        #this number is by pairs to get read num we times 2
        self._rescued["total"] += (new_package.rescued*2)    

        if new_package.results: #This ensures that there are leftover loners
            self._pyramid_idle_counts[loner_type] = 0
            
            try:
                pyramid = self._loner_pyramid[loner_type]
            except KeyError:
                self._loner_pyramid[loner_type] = [[]]
                pyramid = self._loner_pyramid[loner_type]

            if len(pyramid) >= level:
                pyramid.append([])

            for loner_level,loner_path in new_package.results:
                pyramid[loner_level].append(loner_path)

        if self._destroy:
            self._stale_count += 1

    def __handle_primary_task__(self,new_package):
        for match_package in new_package.results:
            self.__handle_match_maker__(match_package)

    def __test_primary_tasks__(self,required_paths=10):
        task_sent = False
        if len(self._primary_store) >= required_paths:
            self.__start_primary_task__()
            task_sent = True

        if task_sent:
            self._primary_store = []
            gc.collect()
            time.sleep(1)
          
    def __start_primary_task__(self):        
        self.__start_job__(self._primary_store,"Primary",-1)

    def __start_matchmaker_task__(self,paths,loner_type,level):
        self.__start_job__(paths,loner_type,level)

    def __start_job__(self,paths,loner_type,level):
        job = JobPackage(paths,loner_type,level)
        self._pending_jobs += 1
        self._chaser_qu.put(job)

    def __pause_monitor__(self):
        if not self._destroy:
            max_jobs = self._max_jobs
            if self._pending_jobs >= max_jobs and not self._file_readers_paused:
                for qu in self._pause_qus:
                    qu.put(1) #pause
                    self.__wait_for_ack__(qu)
                self._file_readers_paused = True

            elif self._pending_jobs < max_jobs and self._file_readers_paused:
                for qu in self._pause_qus:
                    qu.put(0) #unpause
                    self.__wait_for_ack__(qu)
                self._file_readers_paused = False

    def __wait_for_ack__(self,qu):
        count = 0
        while True:
            if count > 10:
                return
            try:
                ack = qu.get(False)
                if ack == 2:
                    return
                else:
                    qu.put(ack)
                    time.sleep(1)
            except Queue2.Empty:
                time.sleep(1)
            count += 1

    def __periodic_action__(self,iterations):

        self.__pause_monitor__()
        chaser_debug = False

        if not self._destroy:
            idle_threshold = 500
            required_paths = 10
        else:
            idle_threshold = 75
            required_paths = 1

        pyramid_idle_counts = self._pyramid_idle_counts 
        empty = True

        for loner_type,pyramid in self._loner_pyramid.items():
            for i,sub_pyramid in enumerate(reversed(pyramid)):
                level = len(pyramid) - (i+1)
                paths = self.__get_paths__(loner_type,sub_pyramid,level)
                if len(paths) > 0:
                    empty = False
                    self.__start_matchmaker_task__(paths,loner_type,level)
                    pyramid_idle_counts[loner_type] = 0
                    for path in paths:#remove sent paths
                        sub_pyramid.remove(path)
                    break
                if len(sub_pyramid) > 0:
                    empty = False
                    pyramid_idle_counts[loner_type] += 1
                del paths

            if pyramid_idle_counts[loner_type] >= idle_threshold:
                idle_success,idle_paths,idle_levels = self.__idle_routine__(pyramid,loner_type)

                #delete after idle success or when processing has finished and idle fails
                if self.__clear_subpyramid__(idle_success, self._destroy,self._pending_jobs):
                    pyramid_idle_counts[loner_type]=0
                    for idle_level,idle_path in izip(idle_levels,idle_paths):
                        if not idle_success:
                            try:
                                self._loner_purgatory[loner_type].append(idle_path)
                            except KeyError:
                                self._loner_purgatory[loner_type] = [idle_path]
                            pyramid[idle_level].remove(idle_path)
                        else:
                            pyramid[idle_level].remove(idle_path)

                del idle_paths,idle_levels

        self.__test_primary_tasks__(required_paths=required_paths)

        if self._destroy:
            self.__post_destroy_report__()
            if not self._rescued["total"] == self._prev_rescued:
                self._stale_count = 0
            self._prev_rescued = self._rescued["total"]
            saved = self.__save_purgatory_loners__()

            if (empty and self._pending_jobs == 0) or self._stale_count > 100:
                finished = True
                if not empty: #essentialy `if stale count`
                    self.__tidy_pyramid__()
                else: 
                    if saved > 0:
                        finished = False
                    if not self.__is_queue_empty__():
                        finished = False
                                       
                if finished:
                    self._finished = True

        if iterations % 20 == 0:
            gc.collect()

        if chaser_debug and iterations % 30 == 0:
            sys.stdout.write("\r %d/%d=%.4f %.2fGB | Proc:%d Des:%d Emp:%d Purg:%d Stale:%d Jobs:%d " %\
                (self._rescued["total"],
                self._total_loners,
                float(self._rescued["total"]+1)/(self._total_loners+1),
                float((self._total_loners-self._rescued["total"]) * 130) / (10**9),
                self._processing,
                self._destroy,
                empty,
                len(self._loner_purgatory),
                self._stale_count,
                self._pending_jobs))
            sys.stdout.flush()

    def __is_queue_empty__(self):
        try:
            pack = self._inqu.get(False,60)
            self._inqu.put(pack)
            return False
        except Queue2.Empty:
            return True

    def __save_purgatory_loners__(self):
        saved = 0
        saved_keys = []
        for loner_type,paths in self._loner_purgatory.items():
            if len(paths) > 1:
                saved += 1
                saved_keys.append(loner_type)
                self._pyramid_idle_counts[loner_type] = 100
                for path in paths:
                    self._loner_pyramid[loner_type][0].append(path)
        for key in saved_keys:
            del self._loner_purgatory[key]
        gc.collect()
        return saved

    def __get_paths__(self,loner_type,sub_pyramid,level):
        if not self._destroy and "UM" in loner_type:
            return []
        task_size = self.__get_task_size__(level,loner_type)
        if len(sub_pyramid) >= task_size:
            if "XX" in loner_type:
                return sub_pyramid[:task_size]
            else:
                return random.sample(sub_pyramid,task_size)
        return []
    
    def __tidy_pyramid__(self):
        for loner_type,pyramid in self._loner_pyramid.items():
            for sub_pyramid in pyramid:
                for x in xrange(len(sub_pyramid)):
                    os.remove(sub_pyramid.pop())
        gc.collect()
            
    def __clear_subpyramid__(self,success,destroy,running):
        if success:
            #If the idle task was sent, clear pyramid
            return True
        elif destroy and not success:
            if running == 0:
                #if main process has finished and no tasks run
                #clear the pyramid
                return True
            else:
                #Tasks still running. Don't clear
                return False 
        else:
            #Idle task not sent, don't clear pyramid
            return False

    def __idle_routine__(self,pyramid,loner_type):
        if "UM" in loner_type and not self._destroy:
            #Don't conduct idle on UMs because we don't expect to see pairs
            #until after processing has finished
            return True,[],[]

        all_paths = []
        for level,sub_pyramid in enumerate(pyramid):
            all_paths.extend(zip([level]*len(sub_pyramid),sub_pyramid))

        success = False
        idle_paths = []
        idle_levels = []
        if len(all_paths) > 1:
            path_n = random.randint(0,len(all_paths))
            if path_n < 2:
                #This gives 3 chances to roll a 2
                #2 paths means a depth heavy search
                path_n = 2

            idle_tuples = random.sample(all_paths,path_n)
            idle_levels,idle_paths = zip( *idle_tuples )

            self.__start_matchmaker_task__(idle_paths,loner_type,max(idle_levels))
            success = True
            del all_paths,idle_tuples
        elif len(all_paths) == 1:
            idle_levels,idle_paths = zip(*all_paths)

        return success,idle_paths,idle_levels

    def __get_task_size__(self,level,loner_type):
        if level % 2 == 0:
            task_size = 2
        else:
            task_size = 3
        return task_size

    def __handler_exit__(self,**kwargs):
        for i in xrange(len(self._chaser_tasks)+1):
           self._chaser_qu.put(DestroyPackage())
        time.sleep(2)
        self._chaser_qu.close()
        del self._chaser_qu
        del self._chaser_tasks
        
        if self._constants.verbose:
            if self._total_loners - self._rescued["total"] == 0:
                self.__standard_output__("\n[Status] All reads succesfully paired") 
            else:
                self.__standard_output__("\n[Status] Couldn't find pairs for %d reads" %\
                    (self._total_loners - self._rescued["total"],))
        for qu in self._pause_qus:
            qu.close()

    def __post_destroy_report__(self):
        if self._constants.verbose:
            if self._post_destroy_count == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()
            if self._post_destroy_count % 5 == 0:
                sys.stdout.write("\r[Status] Read pairing still in progress: %.3f%% complete  " %\
                            ((float(self._rescued["total"]+1) / (self._total_loners+1))*100,))
                sys.stdout.flush()
            self._post_destroy_count += 1

class ChaserTask(Process):
    def __init__(self,object constants,object inqu, object outqu, 
                      int max_reads, dict child_pack):
        super(ChaserTask,self).__init__()
        self._constants = constants
        self._inqu = inqu
        self._outqu = outqu
        self._child_pack = child_pack
        self._max_reads = max_reads

    def run(self):
        primary_handler = PrimaryHandler(constants = self._constants,
                                         max_reads = self._max_reads,
                                         child_pack = self._child_pack,
                                         pid = self.pid)

        match_handler = MatchMakerHandler(constants = self._constants,
                                         max_reads = self._max_reads,
                                         child_pack = self._child_pack,
                                         pid = self.pid)

        while True:
            try:
                package = self._inqu.get(False)
                if type(package) == JobPackage:
                    if package.level == -1:
                        results = primary_handler.run(package.paths)
                    else:
                        results = match_handler.run(package.paths,
                                          package.loner_type,
                                          package.level)
                    self._outqu.put(results)

                elif type(package) == DestroyPackage:
                    break
                
                time.sleep(0.005)

            except Queue2.Empty:
                time.sleep(5)
        return

class ChaserHandler(object):

    def __init__(self,object constants,int max_reads, dict child_pack,int pid):
        self._constants = constants
        self._max_reads = max_reads
        self._child_pack = child_pack
        self._parent_bam = self._child_pack["parent_bam"]
        self._header = self._parent_bam.header
        self._record = 0
        self.pid = pid

    def __get_loner_temp_path__(self,loner_type,level=0,unique=""):
        self._record += 1
        return "%s/%s_%d%d_%s%s.bam" %\
            (self._constants.temp_dir,self.pid,self._record,level,loner_type,unique)

class PrimaryHandler(ChaserHandler):

    def __init__(self,object constants, int max_reads, dict child_pack,int pid):
        super(PrimaryHandler,self).__init__(constants,max_reads,child_pack,pid)
        self.match_handler = MatchMakerHandler(constants,max_reads,child_pack,pid)
        
    def run(self,unsorted_paths):
        #speedup
        write_loner = self.__write_loner__
        #speedup

        loner_holder = self.__get_loner_holder__()
        loner_file_counts = Counter()

        loner_total = 0

        chrom_bins = self._child_pack["chrom_bins"]

        for read in self.__read_generator__(unsorted_paths):
            
            loner_type = self.__get_loner_type__(read,chrom_bins)
            loner_total += 1
            write_loner(loner_holder,loner_file_counts,loner_type,read)

        loner_paths = self.__get_paths_and_close__(loner_holder)

        match_maker_results = []
        for loner_type,paths in loner_paths.items():
            for path in paths:
                match_maker_results.append(self.__start_matchmaker_task__([path],loner_type,-1))

        del loner_paths,loner_file_counts,loner_holder
        return PrimaryResults(match_maker_results,"primary")

    def __start_matchmaker_task__(self,paths,loner_type,level):
        return self.match_handler.run(paths,loner_type,level)
        
    def __read_generator__(self,unsorted_paths):
        for path in unsorted_paths:
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
            os.remove(path)
        gc.collect()

    def __get_reference_id_name__(self,read,bins):
        if read.reference_id == read.next_reference_id:
            return "XX%dv%d" % (read.reference_id,read.next_reference_id)
        else:
            class_bins = map(lambda x : bins[x],sorted((read.reference_id,read.next_reference_id,)))
            #reads on different chromosome
            return "MM%sv%s" % tuple(class_bins)

    def __get_loner_type__(self,read,bins):
        if not read.is_unmapped and not read.mate_is_unmapped:
            return self.__get_reference_id_name__(read,bins)
        elif read.is_unmapped and read.mate_is_unmapped:
            #both unmapped 
            return "UU"
        else:
            #one unmapped read, one mapped read
            return "UM"
            
    def __get_loner_holder__(self):
        holder = {}
        return holder

    def __get_loner_object__(self,loner_type,unique):
        path = self.__get_loner_temp_path__(loner_type,0,unique="-%d" % (unique,) )
        return pysam.AlignmentFile(path,"wb",header=self._header)

    def __write_loner__(self,loner_objects,loner_file_counts,loner_type,read):
        try:
            loner_path = loner_objects[loner_type][-1].filename
            if loner_file_counts[loner_path] == ((self._max_reads / 2) - 1): #Ensure fewer reads than max
                new_bam_object = self.__get_loner_object__(loner_type,len(loner_objects[loner_type]))
                loner_objects[loner_type][-1].close()
                
                loner_objects[loner_type].append(new_bam_object) 
        except KeyError:
            new_bam_object = self.__get_loner_object__(loner_type,0)
            loner_path = new_bam_object.filename
            loner_objects[loner_type] = [new_bam_object]
        
        loner_file_counts[loner_path] += 1
        loner_objects[loner_type][-1].write(read)

    def __get_paths_and_close__(self,loner_holder):   
        paths = {}
        for loner_type,loner_objects in loner_holder.items():
            paths[loner_type] = []
            for loner_object in loner_objects:
                paths[loner_type].append(loner_object.filename)
            loner_object.close()#This close prevents too many open files
        return paths

class MatchMakerHandler(ChaserHandler):

    def __init__(self,object constants,int max_reads, dict child_pack,int pid):
        super(MatchMakerHandler,self).__init__(constants,max_reads,child_pack, pid)
        self.LonerInfo = namedtuple("LonerInfo","paths level loner_type")

    def run(self,loner_paths,loner_type,level):

        loner_info = self.LonerInfo(loner_paths,level,loner_type)
        
        loner_path = self.__get_loner_temp_path__(loner_info.loner_type,loner_info.level+1)
        loner_object = self.__get_loner_object__(loner_path)

        cdef dict pairs = self.__first_pass__(loner_info)
        cdef dict read_pairs = {}
        cdef dict send_pairs = {}
        cdef int rescued = 0
        cdef int written = 0
        cdef list leftover_paths = []

        for read in self.__read_generator__(loner_info,second_pass=True,
                                            leftover_paths=leftover_paths):
            try:
                status = pairs[read.qname]
                try:
                    read_pairs[read.qname].append(read)
                except KeyError:
                    read_pairs[read.qname] = [read]

                if len(read_pairs[read.qname]) == 2:
                    send_pairs[read.qname] = read_pairs[read.qname]
                    del pairs[read.qname]
                    del read_pairs[read.qname]

                #If over hard limit, send to main stream.
                if len(send_pairs) >= 100000:
                    rescued += self.__launch_child_task__(send_pairs)
                    send_pairs = {}
                    gc.collect()

            except KeyError:
                written += 1
                loner_object.write(read)
        loner_object.close()

        if len(send_pairs) > 0:
            rescued += self.__launch_child_task__(send_pairs)

        del read_pairs   
        del pairs

        gc.collect()

        return_paths = []

        if written == 0:
            os.remove(loner_path)
        else:
            return_paths.append((loner_info.level+1,loner_path))

        if len(leftover_paths) > 0:
            for leftover_path in leftover_paths:
                return_paths.append((loner_info.level,leftover_path))

        loner_pack = MatchMakerResults(loner_type=loner_info.loner_type,results=return_paths,
                                       level=loner_info.level+1,chaser_type="match_maker",
                                       rescued=rescued)
        return loner_pack

    def __first_pass__(self,loner_info):
        cdef dict pairs = {}
        for read in self.__read_generator__(loner_info):
            try:
                status = pairs[read.qname]
                pairs[read.qname] = True
            except KeyError:
                pairs[read.qname] = False
        return self.__prepare_for_second_pass__(pairs)

    def __prepare_for_second_pass__(self,pairs):
        return dict(filter( lambda tup : tup[1], iter(pairs.items()) ))

    def __read_generator__(self,loner_info,second_pass=False,leftover_paths=None):
        count = self._max_reads // len(loner_info.paths)
        count_limit = False
        for path in loner_info.paths:
            loner_object = pysam.AlignmentFile(path,"rb")
            bam_iterator = loner_object.fetch(until_eof=True)
           
            for i,read in enumerate(bam_iterator):
                yield read
                if i == count:
                    count_limit = True
                    break 

            if second_pass:
                if count_limit:
                    count_limit = False
                    leftover_count = 0
                    leftover_object = self.__get_leftover_object__(loner_info,len(leftover_paths))
                    for read in bam_iterator:
                        leftover_count += 1
                        leftover_object.write(read)
                    leftover_object.close()
                    
                    if leftover_count > 0:
                        leftover_paths.append(leftover_object.filename)
                    else:
                        os.remove(leftover_object.filename)
                os.remove(path)
            loner_object.close()
        gc.collect()

    def __launch_child_task__(self,pairs):
        child_pack = self._child_pack
        results = child_pack["chaser_task"].handle_pair_dict(pairs,int("%d%d" % (self.pid,self._record,)))
        child_pack["queue"].put(parabam.core.Package(results=results))
        self._record += 1
        return len(pairs)

    def __get_leftover_object__(self,loner_info,unique):
        path = self.__get_loner_temp_path__(loner_info.loner_type,loner_info.level,"-%d" % (unique,))
        leftover = pysam.AlignmentFile(path,"wb",header=self._header)
        return leftover

    def __get_loner_object__(self,path):
        loner_object = pysam.AlignmentFile(path,"wb",header=self._header)
        return loner_object

#.....happily ever after
