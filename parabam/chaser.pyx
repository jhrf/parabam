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
    def __init__(self,paths,loner_type,level,fragment_leftovers=False):
        self.paths = paths
        self.loner_type = loner_type 
        self.level = level
        self.fragment_leftovers = fragment_leftovers

class Handler(parabam.core.Handler):

    def __init__(self,object parent_bam, object output_paths,
                      object inqu,object constants,object pause_qus,
                      dict out_qu_dict,object Task):

        super(Handler,self).__init__(parent_bam = parent_bam,
                                     output_paths = output_paths,
                                     inqu=inqu,
                                     constants=constants,
                                     pause_qus=pause_qus,
                                     out_qu_dict=out_qu_dict,report = False)
        
        class RuleTask(Task):
            def __init__(self,parent_bam,constants):
                super(RuleTask,self).__init__(parent_bam=parent_bam,
                                                    inqu=None,
                                                    outqu=None,
                                                    statusqu=None,
                                                    task_size=0,
                                                    constants=constants)

            def handle_pair_dict(self,pairs,unique):
                self.unique = unique
                results = self.__generate_results__(pairs)
                results["Total"] = 0
                self._dealt += 1
                time.sleep(.01)
                return results

            def __process_task_set__(self,iterator):
                rule = self._rule
                parent_bam = self._parent_bam
                handle_output = self.__handle_rule_output__
                user_constants = self._constants.constants

                for read_name,pair in iterator.items():
                    rule_output = rule(pair,user_constants,parent_bam)
                    handle_output(rule_output,pair)

            def __get_temp_path__(self,identity):
                file_name = "%s_%d_%d_%s" %\
                    (identity,self.unique,self._dealt, os.path.split(\
                                                self._parent_bam.filename)[1])
                return os.path.join(self._temp_dir,file_name)

        self._loner_pyramid = self.__instalise_loner_pyramid__()
        self._pyramid_idle_counts = Counter()
        self._loner_purgatory = {}
        self._files_in_pyramid = 0

        self._periodic_interval = 5
        self._total_loners = 0
        self._rescued = Counter()
        self._prev_rescued = 0
        self._stale_count = 0
        self._stale_thresh = 250

        self._child_pack = {"parent_bam":parent_bam,
                            "queue":out_qu_dict["main"],
                            "chaser_task":RuleTask(parent_bam,
                                                     self._constants)}

        self._chaser_qu = Queue()
        self._chaser_tasks = self.__create_chaser_tasks__(\
                                                self._constants.total_procs)
        self._pending_jobs = 0
        self._max_jobs = (self._constants.total_procs * 3)
        self._file_readers_paused = False

        self._post_destroy_count = 0
        self._post_destroy_thresh = 50

        #self._print_chaser_debug = (lambda iterations: False)
        self._print_chaser_debug = self.__print_debug__

    def __print_debug__(self,iterations):
        if iterations % 50 == 0 or iterations == -1:

            if iterations == -1:
                sys.stdout.write("\nFINISHED\n")

            sys.stdout.write(("\r %d/%d=%.4f %.2fGB | Proc:%d Des:%d "
                         "Purg:%d Stale(Thresh):%d(%d) Jobs:%d Files:%d  ") %\
                (self._rescued["total"],
                self._total_loners,
                float(self._rescued["total"]+1)/(self._total_loners+1),
                float((self._total_loners\
                            - self._rescued["total"]) * 130) / (10**9),
                self._processing,
                self._destroy,
                len(self._loner_purgatory),
                self._stale_count,
                self._stale_thresh,
                self._pending_jobs,
                self._files_in_pyramid))

            if iterations == -1:
                sys.stdout.write("\nFINISHED\n")
            sys.stdout.flush()

        if iterations % 2500 == 0:
            print ""

    def __create_chaser_tasks__(self,total_tasks):
        tasks = []
        for i in xrange(total_tasks):
            task = ChaserTask(constants=self._constants,
                              inqu=self._chaser_qu,
                              outqu=self._inqu,
                              max_reads=4000000,
                              child_pack=self._child_pack)
            task.start()
            tasks.append(task)
        return tasks
           
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
        if chaser_type == "match_maker":
            self.__handle_match_maker__(new_package)
            self._pending_jobs -= 1
        elif chaser_type == "origin":
            self.__handle_origin_task__(new_package)
        else:
            print "Error!"

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
                self.__add_to_pyramid__(loner_type,loner_level,loner_path)

        if self._destroy: 
            #TODO:should this move to new_package_action?
            #     think carefully before doing this as this will
            #     change the way that stale counting behaves
            self._stale_count += 1

    def __handle_origin_task__(self,new_package):
        for loner_count, primary_paths  in new_package.results:
            if loner_count > 0:
                self._total_loners += loner_count
                for loner_type,loner_path in primary_paths.items():

                    match_package = MatchMakerResults(
                                          loner_type=loner_type,
                                          results=((0,loner_path),),
                                          chaser_type="match_maker",
                                          rescued=0,
                                          level=0)
                    self.__handle_match_maker__(match_package)

    def __start_matchmaker_task__(self,paths,loner_type,level):
        self.__start_job__(paths,loner_type,level)

    def __start_job__(self,paths,loner_type,level):
        fragment_leftovers = self.__request_fragment_output__(paths)

        job = JobPackage(paths,loner_type,level,fragment_leftovers)
        self._pending_jobs += 1
        self._chaser_qu.put(job)

    def __request_fragment_output__(self,paths):
        if self._destroy and self._stale_count > 10:
            has_large_file = any([ (os.path.getsize(path) / 100) \
                                             > 1000000 for path in paths ])
            if has_large_file:
                if self._stale_thresh < 250:
                    self._stale_thresh = 250
                return True
        return False

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

        if not self._destroy:
            idle_threshold = 500
        else:
            idle_threshold = 200

        pyramid_idle_counts = self._pyramid_idle_counts 

        for loner_type,pyramid in self._loner_pyramid.items():
            for i,sub_pyramid in enumerate(reversed(pyramid)):
                level = len(pyramid) - (i+1)
                paths = self.__get_paths__(loner_type,sub_pyramid,level)
                if len(paths) > 0:
                    self.__start_matchmaker_task__(paths,loner_type,level)
                    pyramid_idle_counts[loner_type] = 0
                    for path in paths:#remove sent paths
                        self.__remove_from_pyramid__(loner_type,level,path)
                    break
                if len(sub_pyramid) > 0:
                    pyramid_idle_counts[loner_type] += 1
                del paths

            if pyramid_idle_counts[loner_type] >= idle_threshold:
                self.__test_pyramid_idle_count__(pyramid,
                                                  loner_type,
                                                  pyramid_idle_counts)

        self._print_chaser_debug(iterations)
        self.__finish_test__()

        if iterations % 100 == 0:
            gc.collect()

    def __test_pyramid_idle_count__(self,pyramid,loner_type,pyramid_idle_counts):
        #Idle counting is the mechanism whereby loners that have no
        #match can be entered into `purgatory` thus clearing the pyramid.
        #This mechanism catches the case where there are genuinley lonerred
        #reads in the BAM file.

        idle_success,idle_paths,idle_levels = \
                            self.__idle_routine__(pyramid,loner_type)

        # Delete after idle success or when 
        # processing has finished and idle fails
        if self.__clear_subpyramid__(idle_success, 
                                     self._destroy,
                                     self._pending_jobs):
            
            pyramid_idle_counts[loner_type]=0
            for idle_level,idle_path in izip(idle_levels,idle_paths):
                if not idle_success:
                    try:
                        self._loner_purgatory[loner_type].\
                                                append(idle_path)
                    except KeyError:
                        self._loner_purgatory[loner_type] = [idle_path]
                self.__remove_from_pyramid__(loner_type,idle_level,idle_path)
                    
        del idle_paths,idle_levels

    def __finish_test__(self):
        if self._destroy:
            self.__post_destroy_report__()
            if not self._rescued["total"] == self._prev_rescued:
                self._stale_count = 0
                self.__set_new_stale_thresh__()

            self._prev_rescued = self._rescued["total"]
            saved = self.__save_purgatory_loners__()

            finished = False
            if self.__pyramid_is_empty__(saved):
                finished = True
            elif self._stale_count > self._stale_thresh:
                finished = self.__tidy_pyramid__()

            if finished:
                self._print_chaser_debug(-1)
                self._finished = True

    def __pyramid_is_empty__(self,saved):
        return (self._files_in_pyramid == 0 and \
                self._pending_jobs == 0 and \
                saved == 0 and \
                self.__is_queue_empty__())

    def __set_new_stale_thresh__(self):
        self._stale_thresh = (10 + (self._files_in_pyramid * 5) \
                                 + (self._pending_jobs * 8))

    def __is_queue_empty__(self):
        try:
            pack = self._inqu.get(True,10)
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
                    self.__add_to_pyramid__(loner_type,0,path)
        for key in saved_keys:
            del self._loner_purgatory[key]
        gc.collect()
        return saved

    def __add_to_pyramid__(self,loner_type,level,path):
        self._loner_pyramid[loner_type][level].append(path)
        self._files_in_pyramid += 1

    def __remove_from_pyramid__(self,loner_type,level,path):
        self._loner_pyramid[loner_type][level].remove(path)
        self._files_in_pyramid -= 1

    def __get_paths__(self,loner_type,sub_pyramid,level):
        if not self._destroy and "UM" in loner_type:
            return []
        task_size = self.__get_task_size__(level,loner_type)
        if len(sub_pyramid) >= task_size:
            if "XX" in loner_type and not self._destroy:
                return sub_pyramid[:task_size] 
            else:
                return random.sample(sub_pyramid,task_size)
        return []
    
    def __tidy_pyramid__(self):
        reads_found = self.__wait_for_remaining_jobs__()

        if not reads_found:
            for loner_type,pyramid in self._loner_pyramid.items():
                for level,sub_pyramid in enumerate(pyramid):
                    for path in sub_pyramid:
                        self.__remove_from_pyramid__(loner_type,level,path)
            gc.collect()
            return True
        else:
            return False

    def __wait_for_remaining_jobs__(self):
        def waiting_update(found,total):
            sys.stdout.write(("\r\t- Waiting for tasks to finish: %d/%d ") % \
                                                            (found,total))
            sys.stdout.flush()

        if self._constants.verbose:
            sys.stdout.write("\n")
            sys.stdout.flush()

        prev_found = self._rescued["total"]
        jobs_at_start = self._pending_jobs
        jobs_found = 0

        i = 0
        while self._pending_jobs > 0:
            if self._constants.verbose and i % 50 == 0:
                jobs_at_start = self._pending_jobs
                waiting_update(jobs_found,jobs_at_start)
            try:
                pack = self._inqu.get(True,10)
                self.__new_package_action__(pack)
                jobs_found += 1
            except Queue2.Empty:
                time.sleep(1)
            i += 1

        if self._constants.verbose:
            waiting_update(jobs_found,jobs_at_start)
            sys.stdout.write("\n")
            sys.stdout.flush()

        print "Saved from wait %d" % (self._rescued["total"] - prev_found)
        return prev_found < self._rescued["total"]

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
        if level == 0:
            task_size = 10
        elif level % 2 == 0:
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
                sys.stdout.write(\
                     "\r\t- Read pairing in progress: 100.00% complete\n")
                sys.stdout.flush()
            self.__standard_output__("\t- Unpaired reads: %d" %\
                (self._total_loners - self._rescued["total"],))
        for qu in self._pause_qus:
            qu.close()

    def __post_destroy_report__(self):
        if self._constants.verbose:
            if self._post_destroy_count >= self._post_destroy_thresh:
                if self._post_destroy_count % 100 == 0 and\
                   self._constants.verbose == 2:
                    sys.stdout.write(\
                        "\r\t- Read pairing in progress: %.2f%% complete  " %\
                          ((float(self._rescued["total"]+1) / (self._total_loners+1))*100,))
                    sys.stdout.flush()
                elif self._post_destroy_count % 1000 == 0 and\
                     self._constants.verbose ==1:
                    sys.stdout.write(\
                        "\n\t- Read pairing in progress: %.2f%% complete" %\
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
        match_handler = MatchMakerHandler(constants = self._constants,
                                         max_reads = self._max_reads,
                                         child_pack = self._child_pack,
                                         pid = self.pid)

        while True:
            try:
                package = self._inqu.get(False)
                if type(package) == JobPackage:
                    results = match_handler.run(package.paths,
                                      package.loner_type,
                                      package.level,
                                      package.fragment_leftovers)
                    self._outqu.put(results)

                elif type(package) == DestroyPackage:
                    break
                
                time.sleep(0.005)

            except Queue2.Empty:
                time.sleep(5)
        return

class MatchMakerHandler(object):

    def __init__(self,
                 object constants,
                 int max_reads,
                 dict child_pack,
                 int pid):

        self._constants = constants
        self._max_reads = max_reads
        self._child_pack = child_pack
        self._parent_bam = self._child_pack["parent_bam"]
        self._header = self._parent_bam.header
        self._record = 0
        self.pid = pid

        self.LonerInfo = namedtuple("LonerInfo","paths level loner_type")

    def run(self,loner_paths,loner_type,level,fragment_leftovers):

        #TODO: This is pretty old code... Not sure if the memory mimimisation
        #      of the two pass system is worth it...
        #      Consider picking through and simplifying here

        loner_info = self.LonerInfo(loner_paths,level,loner_type)
        
        loner_path = self.__get_loner_temp_path__(\
                                    loner_info.loner_type,loner_info.level+1)
        loner_object = self.__get_loner_object__(loner_path)

        cdef dict pairs = self.__first_pass__(loner_info,fragment_leftovers)
        cdef dict read_pairs = {}
        cdef dict send_pairs = {}
        cdef int rescued = 0
        cdef int written = 0
        cdef list leftover_paths = []

        for read in self.__read_generator__(loner_info,second_pass=True,
                                            leftover_paths=leftover_paths,
                                            fragment_leftovers=fragment_leftovers):
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
 
        loner_pack = MatchMakerResults(loner_type=loner_info.loner_type,
                                       results=return_paths,
                                       level=loner_info.level+1,
                                       chaser_type="match_maker",
                                       rescued=rescued)
        return loner_pack

    def __get_loner_temp_path__(self,loner_type,level=0,unique=""):
        self._record += 1
        return "%s/%s_%d_%d_%s_%s.bam" %\
            (self._constants.temp_dir,self.pid,
             self._record,level,loner_type,unique)

    def __first_pass__(self,loner_info,fragment_leftovers):
        cdef dict pairs = {}
        for read in self.__read_generator__(loner_info,
                                      fragment_leftovers=fragment_leftovers):
            try:
                status = pairs[read.qname]
                pairs[read.qname] = True
            except KeyError:
                pairs[read.qname] = False
        return self.__prepare_for_second_pass__(pairs)

    def __prepare_for_second_pass__(self,pairs):
        return \
            dict( (qname, status,) for (qname,status) in pairs.items() if status)

    def __read_generator__(self,loner_info,
                            second_pass=False,leftover_paths=None,
                            fragment_leftovers=False):

        if not fragment_leftovers:
            count = self._max_reads // len(loner_info.paths)
        else:
            count = 50000

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
                    self.__write_leftovers__(bam_iterator,
                                             fragment_leftovers,
                                             loner_info,
                                             leftover_paths)

                os.remove(path)
            loner_object.close()
        gc.collect()

    def __write_leftovers__(self,bam_iterator,fragment_leftovers,
                                 loner_info,leftover_paths):

        leftover_file = self.__get_leftover_object__(\
                                        loner_info,len(leftover_paths))
        leftover_count = 0
        for read in bam_iterator:
            if fragment_leftovers and\
               leftover_count % 25000 == 0 and\
               leftover_count > 0:
                
                self.__add_to_leftover_paths__(leftover_file,
                                               leftover_paths)

                leftover_file = self.__get_leftover_object__(\
                                        loner_info,len(leftover_paths))
 
            leftover_file.write(read)
            leftover_count += 1

        if leftover_count > 0:
            self.__add_to_leftover_paths__(leftover_file,
                                           leftover_paths)
        else:
            empty_leftover_path = leftover_file.filename
            leftover_file.close()
            os.remove(empty_leftover_path)

    def __add_to_leftover_paths__(self,leftover_file,leftover_paths):
        leftover_paths.append(leftover_file.filename)
        leftover_file.close()

    def __launch_child_task__(self,pairs):
        child_pack = self._child_pack
        results = child_pack["chaser_task"].handle_pair_dict(pairs,\
                                    int("%d%d" % (self.pid,self._record,)))
        child_pack["queue"].put(parabam.core.Package(results=results))
        self._record += 1
        return len(pairs)

    def __get_leftover_object__(self,loner_info,unique):
        path = self.__get_loner_temp_path__(loner_info.loner_type,\
                                        loner_info.level,"lftovr%d" % (unique,))
        leftover = pysam.AlignmentFile(path,"wb",header=self._header)
        return leftover

    def __get_loner_object__(self,path):
        loner_object = pysam.AlignmentFile(path,"wb",header=self._header)
        return loner_object

#.....happily ever after
