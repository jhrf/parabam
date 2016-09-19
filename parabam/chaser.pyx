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
from pprint import pprint as ppr

class LonerRecord(object):
    def __init__(self, 
                 loner_type, 
                 source, 
                 target, 
                 level, 
                 temp_dir,
                 salt,
                 header):
        
        self.loner_type = loner_type
        self.source = source
        self.target = target
        self.level = level
        self.temp_dir = temp_dir
        self.salt = salt
        self.header = header

        self.count = 0
        self.path = self.__get_path__()
        self.file_object = None

    def __get_path__(self):
        basename = "loner_%s_%s_%s-%d.bam" % (self.loner_type,
                                            self.source,
                                            self.salt,
                                            self.level)
        path = os.path.join(self.temp_dir,basename)
        return path

    def file_exists(self):
        return os.path.exists(self.path)

    def new_record(self, salt, level_mod):
        return LonerRecord(
            loner_type = self.loner_type,
            source = self.source,
            target = self.target,
            level = self.level + level_mod,
            temp_dir = self.temp_dir,
            salt = salt,
            header = self.header)

    def write(self, read):
        if self.file_object is None:
            if os.path.exists(self.path):
                print 'WARNING: Attempt to overwrite loner file! Data lost!'
            self.file_object = pysam.AlignmentFile(self.path,
                                                   "wb",
                                                   header=self.header)
        self.file_object.write(read)
        self.count += 1

    def read(self):
        if self.file_object is None:
            self.file_object = pysam.AlignmentFile(self.path,"rb")
            self.bam_iter = self.file_object.fetch(until_eof=True)
        return self.bam_iter.next()

    def close(self):
        if self.file_object is not None:
            self.file_object.close()

        self.file_object = None
        self.bam_iter = None

class ChaserResults(parabam.core.Package):
    def __init__(self,results,chaser_type):
        super(ChaserResults,self).__init__(results)
        self.chaser_type=chaser_type

class MatchMakerResults(ChaserResults):
    def __init__(self,loner_type,results,chaser_type,rescued,level):
        super(MatchMakerResults,self).__init__(results,chaser_type)
        self.rescued = rescued
        self.loner_type = loner_type
        self.level = level

class JobPackage(object):
    def __init__(self,records,level,fragment_leftovers=False):
        self.records = records
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
        
            try:
                pyramid = self._loner_pyramid[loner_type]
            except KeyError:
                self._loner_pyramid[loner_type] = [[]]
                pyramid = self._loner_pyramid[loner_type]

            if len(pyramid) >= level:
                pyramid.append([])

            for loner_record in new_package.results:
                self.__add_to_pyramid__(loner_record)

        if self._destroy: 
            #TODO: should this move to new_package_action?
            #      think carefully before doing this as this will
            #      change the way that stale counting behaves
            self._stale_count += 1

    def __handle_origin_task__(self,new_package):
        for primary_paths  in new_package.results:
            for loner_info,record in primary_paths.items():
                self._total_loners += record.count

                match_package = MatchMakerResults(
                                      loner_type=record.loner_type,
                                      results=[record],
                                      chaser_type="match_maker",
                                      rescued=0,
                                      level=0)
                self.__handle_match_maker__(match_package)

    def __start_matchmaker_task__(self, records, loner_type, level):
        self.__start_job__(records, level)
        self._pyramid_idle_counts[loner_type] = 0

    def __start_job__(self, records, level):
        fragment_leftovers = self.__request_fragment_output__(records)

        job = JobPackage(records, level, fragment_leftovers)
        self._pending_jobs += 1
        self._chaser_qu.put(job)

    def __request_fragment_output__(self,records):
        if self._destroy:
            has_large_file = any([ record.count > 1000000 \
                                    for record in records ])
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
                records = \
                    self.__get_records_for_search__(loner_type,sub_pyramid,level)
                if len(records) > 0:
                    self.__start_matchmaker_task__(records,loner_type,level)

                    for record in records:#remove sent records
                        self.__remove_from_pyramid__(record)
                    break

                if len(sub_pyramid) > 0:
                    pyramid_idle_counts[loner_type] += 1
                del records

            if pyramid_idle_counts[loner_type] >= idle_threshold:
                self.__test_pyramid_idle_count__(pyramid,
                                                 loner_type)

        self._print_chaser_debug(iterations)
        self.__finish_test__()

        if iterations % 100 == 0:
            gc.collect()

    def __test_pyramid_idle_count__(self, pyramid, loner_type):
        #Idle counting is the mechanism whereby loners that have no
        #match can be entered into `purgatory` thus clearing the pyramid.
        #This mechanism catches the case where there are genuinley lonerred
        #reads in the BAM file.

        idle_success,idle_records = \
                            self.__idle_routine__(pyramid,loner_type)

        # Delete after idle success or when 
        # processing has finished and idle fails
        if self.__clear_subpyramid__(idle_success, 
                                     self._destroy,
                                     self._pending_jobs):

            for idle_record in idle_records:
                self.__remove_from_pyramid__(idle_record)
                if not idle_success:
                    self.__add_to_purgatory__(idle_record)
                    #Level is zero for all purgatory loners

        del idle_records

    def __add_to_purgatory__(self, record):
        self.__prepare_purgatory__(record)
        record.level = 0
        self._loner_purgatory[record.loner_type][record.source].append(record)

    def __prepare_purgatory__(self, record):

        loner_type = record.loner_type
        if loner_type not in self._loner_purgatory:
            self._loner_purgatory[loner_type] = {}
            self._loner_purgatory[loner_type][record.source] = []
            if not record.source == record.target:
                self._loner_purgatory[loner_type][record.target] = []

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
        if "ccM-U" in loner_type and not self._destroy:
            #Don't conduct idle on UMs because we don't expect to see pairs
            #until after processing has finished
            return True,[]

        all_records = []
        for level,sub_pyramid in enumerate(pyramid):
            all_records.extend(sub_pyramid)

        success = False
        idle_records = []

        if len(all_records) > 1:
        
            if all_records[0].source == all_records[0].target:
                success = True
                path_n = random.randint(0,len(all_records))
                if path_n < 2:
                    #This gives 3 chances to roll a 2
                    #2 paths means a depth heavy search
                    path_n = 2

                idle_records = random.sample(all_records,path_n)
            
            else:
                counts = Counter({all_records[0].source:0,
                                  all_records[1].source:0})
                for record in all_records:
                    counts[record.source] += 1

                success = 0 not in counts.values()
                idle_records = all_records

            if success:
                self.__start_matchmaker_task__(idle_records,
                                               loner_type,
                                               0)
            del all_records

        elif len(all_records) == 1:
            idle_levels = all_records

        return success,idle_records

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
        for loner_type,source_dict in self._loner_purgatory.items():
            counts = [len(records) for source,records in source_dict.items()]
            if 0 not in counts:
                saved += sum(counts)
                saved_keys.append(loner_type)
                self._pyramid_idle_counts[loner_type] = 100
                for source,records in source_dict.items():
                    for record in records:
                        self.__add_to_pyramid__(record)
        for key in saved_keys:
            del self._loner_purgatory[key]
        gc.collect()
        return saved

    def __add_to_pyramid__(self,record):
        self._loner_pyramid[record.loner_type][record.level].append(record)
        self._files_in_pyramid += 1

    def __remove_from_pyramid__(self,record):
        self._loner_pyramid[record.loner_type][record.level].remove(record)
        self._files_in_pyramid -= 1

    def __get_records_for_search__(self,loner_type,sub_pyramid,level):
        records = []

        task_size = self.__get_task_size__(level,loner_type)
        if len(sub_pyramid) >= task_size:
            if sub_pyramid[0].source == sub_pyramid[0].target:
                records =  sub_pyramid[:task_size]
            elif not sub_pyramid[0].source == sub_pyramid[0].target:
                records = self.__get_dif_type_records__(sub_pyramid, task_size)
        return records

    def __get_dif_type_records__(self, sub_pyramid, task_size):

        task_records = []
        records_by_type, type_counts, loner_types =\
                 self.__sperate_record_by_type__(sub_pyramid)      

        missing_type = 0 in type_counts.values()

        if not missing_type:
            if len(sub_pyramid) == task_size:
                task_records = list(sub_pyramid)
            else:
                hi_type,lo_type = zip(*type_counts.most_common(2))[0]

                upper_sample_limit = min(task_size - 1, type_counts[lo_type])

                lo_type_sample_n = random.randint(1,upper_sample_limit)
                print len(sub_pyramid), type_counts[hi_type], type_counts[lo_type]

                lo_sample = random.sample(records_by_type[lo_type],
                                             lo_type_sample_n)
                hi_sample = random.sample(records_by_type[hi_type],
                                             task_size - lo_type_sample_n)

                task_records.extend(lo_sample)
                task_records.extend(hi_sample)
                del lo_sample, hi_sample

        if missing_type and not self._destory:
            self._pyramid_idle_counts[sub_pyramid[0].loner_type] = 0

        return task_records

    def __sperate_record_by_type__(self, sub_pyramid):
        type_1 = sub_pyramid[0].source
        type_2 = sub_pyramid[0].target

        records = {type_1:[],type_2:[]}
        counts = Counter({type_1:0,type_2:0})

        for record in sub_pyramid:
            if record.source == type_1:
                records[type_1].append(record)
                counts[type_1] += 1
            else:
                records[type_2].append(record)
                counts[type_2] += 1

        return records, counts, (type_1, type_2)

    def __tidy_pyramid__(self):
        reads_found = self.__wait_for_remaining_jobs__()

        if not reads_found:
            for loner_type,pyramid in self._loner_pyramid.items():
                for level,sub_pyramid in enumerate(pyramid):
                    for record in sub_pyramid:
                        self.__remove_from_pyramid__(record)
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

    def __get_task_size__(self,level,loner_type):
        if level == 0:
            task_size = 10
        elif level % 2 == 0:
            task_size = 4
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
                     "\r\t- Read pairing in progress: 100.00%% complete\n")
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
                    results = match_handler.run(package.records,
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
        self._unique_tracker = 0
        self.pid = pid

    def run(self, loner_records, level, fragment_leftovers):

        #TODO: This is pretty old code... Not sure if the memory mimimisation
        #      of the two pass system is worth it...
        #      Consider picking through and simplifying here

        unpaired_records = self.__get_unpaired_records__(loner_records,
                                                            level_mod=1)
        
        cdef dict pairs = self.__first_pass__(loner_records,fragment_leftovers)
        cdef dict read_pairs = {}
        cdef dict send_pairs = {}
        cdef int rescued = 0
        cdef list leftover_records = []

        for read, source in self.__read_generator__(
                                        loner_records,
                                        second_pass=True,
                                        leftover_records=leftover_records,
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
                unpaired_records[source].write(read)

        if len(send_pairs) > 0:
            rescued += self.__launch_child_task__(send_pairs)

        del read_pairs   
        del pairs

        return_paths = self.__get_return_paths__\
                            (unpaired_records, leftover_records)

        gc.collect()
 
        loner_pack = MatchMakerResults(loner_type=loner_records[0].loner_type,
                                       results=return_paths,
                                       level=level+1,
                                       chaser_type="match_maker",
                                       rescued=rescued)
        return loner_pack

    def __get_return_paths__(self, loner_records, leftover_records):
        return_paths = []

        for source, record in loner_records.items():
            if record.count > 0:
                return_paths.append(record)
            record.close()

        for record in leftover_records:
            return_paths.append(record)

        return return_paths

    def __get_unpaired_records__(self,records, level_mod):
        unpaired_records = {}

        for record in records:
            if record.source not in unpaired_records:
                new_record = self.__get_new_record__(record, level_mod)
                unpaired_records[record.source] = new_record
        return unpaired_records

    def __get_new_record__(self,model_record,level_mod):
        self._unique_tracker += 1
        salt = "%s-%d" % (self.pid, self._unique_tracker,) 
        return model_record.new_record(salt = salt, level_mod = level_mod)

    def __first_pass__(self,records, fragment_leftovers):
        cdef dict pairs = {}
        for read, source in self.__read_generator__(records,
                                    fragment_leftovers=fragment_leftovers):
            if read.qname in pairs:
                pairs[read.qname] = True
            else:
                pairs[read.qname] = False
        return self.__prepare_for_second_pass__(pairs)

    def __prepare_for_second_pass__(self,pairs):
        return \
            dict( (qname, status,) for (qname,status) in pairs.items() if status)

    def __read_generator__(self,records,
                            second_pass=False,
                            leftover_records=None,
                            fragment_leftovers=False):

        if not fragment_leftovers:
            count = self._max_reads // len(records)
        else:
            count = 50000

        for record in records:

            for i in range(count+1):
                try:
                    read = record.read()
                    yield read, record.source
                except StopIteration:
                    pass

            if second_pass and i == count:
                self.__write_leftovers__(record,
                                         fragment_leftovers,
                                         leftover_records)
            record.close()
            if second_pass:
                os.remove(record.path)

        gc.collect()

    def __write_leftovers__(self,record,
                                fragment_leftovers,
                                leftover_records):

        leftover = self.__get_new_record__(record, level_mod=0)
        
        for read in record.bam_iter:
            if fragment_leftovers and\
               leftover.count % 25000 == 0 and\
               leftover.count > 0:
                
                self.__add_to_leftover_records__(leftover,
                                                 leftover_records)
                leftover = self.__get_new_record__(record, level_mod=0)
            leftover.write(read)

        if leftover.count > 0:
            self.__add_to_leftover_records__(leftover,leftover_records)

    def __add_to_leftover_records__(self, record, leftover_records):
        leftover_records.append(record)
        record.close()

    def __launch_child_task__(self,pairs):
        child_pack = self._child_pack
        results = child_pack["chaser_task"].handle_pair_dict(pairs,\
                                    int("%d%d" % (self.pid,self._unique_tracker,)))
        child_pack["queue"].put(parabam.core.Package(results=results))
        self._unique_tracker += 1
        return len(pairs)

#.....happily ever after
