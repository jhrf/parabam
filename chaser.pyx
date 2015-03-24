#Once upon a time
import pysam
import time
import sys
import os
import gc
import parabam

import Queue as Queue2
import numpy as np

from itertools import izip
from collections import Counter
from multiprocessing import Queue,Process

class ChaserPackage(parabam.core.Package):
    def __init__(self,results,chaser_type):
        super(ChaserPackage,self).__init__(results)
        self.chaser_type=chaser_type

class OriginPackage(ChaserPackage):
    def __init__(self,results,chaser_type):
        super(OriginPackage,self).__init__(results,chaser_type)

class MatchMakerPackage(ChaserPackage):
    def __init__(self,loner_type,results,chaser_type,rescued,level):
        super(MatchMakerPackage,self).__init__(results,chaser_type)
        self.rescued = rescued
        self.loner_type = loner_type
        self.level = level

class Handler(parabam.core.Handler):

    def __init__(self,object parent_bam, object output_paths,
                      object inqu,object constants,object pause_qus,
                      dict out_qu_dict,object Task):

        super(Handler,self).__init__(parent_bam = parent_bam,output_paths = output_paths,
                             inqu=inqu,constants=constants,pause_qus=pause_qus,
                             out_qu_dict=out_qu_dict,report = False)
        
        class ChaserTask(Task):
            def __init__(self,parent_bam,constants):
                super(ChaserTask,self).__init__(parent_bam=parent_bam,
                                                        inqu=None,
                                                        outqu=None,
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

                #StopIteration caught in parabam.core.Task.run
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
                            "chaser_task":ChaserTask(parent_bam,self._constants)}

        self._chaser_task_max = 8 #TODO: work this out properly
        self._tasks = []
        self._primary_store = self.__instalise_primary_store__()

        self._primary_complete = True

        self._max_reads = 5000000


    #Need to refactor here to stop code duplication
    def __wait_for_tasks__(self,list active_tasks,int max_tasks):
        update_tasks = self.__update_tasks__ #optimising alias
        update_tasks(active_tasks)
        cdef int currently_active = len(active_tasks)
        self._active_count = currently_active

        if max_tasks > currently_active:
            return

        if not self._destroy:
            for qu in self._pause_qus:
                qu.put(True)

        while(max_tasks <= currently_active and currently_active > 0):
            update_tasks(active_tasks)
            currently_active = len(active_tasks)
            time.sleep(1)

        if not self._destroy:
            for qu in self._pause_qus:
                qu.put(False)
        return

    def __instalise_primary_store__(self):
        paths = []
        return paths

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

    def __instalise_loner_pyramid__(self):
        loner_pyramid = {}
        return loner_pyramid

    def __new_package_action__(self,new_package,**kwargs):        
        chaser_type = new_package.chaser_type
        if chaser_type == "origin":
            self.__handle_origin_task__(new_package)
        elif chaser_type == "match_maker":
            self.__handle_match_maker__(new_package)
        else:
            print "Error!"
        self.__update_tasks__(self._tasks)

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
        self.__wait_for_tasks__(self._tasks,self._chaser_task_max) #wait for existing tasks
        paths = self._primary_store
        new_task = PrimaryTask(paths,self._inqu,self._constants,self._max_reads,self._child_pack)
        new_task.start()
        self._tasks.append(new_task)

    def __start_matchmaker_task__(self,paths,loner_type,level):
        self.__wait_for_tasks__(self._tasks,self._chaser_task_max)

        new_task = MatchMakerTask(paths,self._inqu,self._constants,
                    loner_type,level,self._max_reads,self._child_pack)
        new_task.start()

        self._tasks.append(new_task)

    def __periodic_action__(self,iterations):
        running = len(self._tasks)

        if not self._destroy:
            idle_threshold = 75
            required_paths = 10
        else:
            idle_threshold = 5
            required_paths=1

        pyramid_idle_counts = self._pyramid_idle_counts 

        empty = True

        for loner_type,pyramid in self._loner_pyramid.items():
                for i,sub_pyramid in enumerate(reversed(pyramid)):
                    level = len(pyramid) - (i+1)
                    paths = self.__get_paths__(loner_type,sub_pyramid,level)
                    if len(paths) > 0:
                        empty = False
                        self.__start_matchmaker_task__(paths,loner_type,level)
                        running += 1
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
                    if self.__clear_subpyramid__(idle_success, self._destroy ,running):
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
                                running += 1

                    del idle_paths,idle_levels

        self.__test_primary_tasks__(required_paths=required_paths)

        if self._destroy:
            if self._rescued["total"] == self._prev_rescued:
                self._stale_count += 1
            else:
                self._stale_count = 0
            self._prev_rescued = self._rescued["total"]
            saved = self.__save_purgatory_loners__()

            if (empty and running == 0) or self._stale_count == 200:
                finished = True
                if not empty:
                    self.__wait_for_tasks__(self._tasks,max_tasks=0)
                    self.__tidy_pyramid__()
                else: 
                    if saved > 0:
                        finished = False
                    if not self.__is_queue_empty__():
                        finished = False
                                       
                if finished:
                    self._finished = True

        gc.collect()

#        if iterations % 30 == 0:
#            sys.stdout.write("\r %d/%d=%.5f | Processing:%d Empty:%d Purgatory:%d Stale:%d Tasks:%d "  %\
#                (self._rescued["total"],
#                self._total_loners,
#                float(self._rescued["total"]+1)/(self._total_loners+1),
#                self._processing,
#                empty,
#                len(self._loner_purgatory),
#                self._stale_count,
#                len(self._tasks)))
#            sys.stdout.flush()

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
        task_size = self.__get_task_size__(level,loner_type)
        paths = []
        if len(sub_pyramid) >= task_size:
            for path in sub_pyramid:
                paths.append(path)
                if len(paths) == task_size:
                    break
        return paths
    
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
        idle_paths = []
        levels = []
        task_size = 5
        for i,sub_pyramid in enumerate(pyramid):
            for path in sub_pyramid:
                idle_paths.append(path)
                levels.append(i)
                if len(idle_paths) == task_size:
                    break 
        success = False
        if len(idle_paths) > 1:
            self.__start_matchmaker_task__(idle_paths,loner_type,levels[-1])
            success = True
        return success,idle_paths,levels
       
    def __get_task_size__(self,level,loner_type):
        if "XX" in loner_type:
            if level % 2 == 0:
                task_size = 2
            else:
                task_size = 3
        elif loner_type == "U-M" and level == 0 and not self._destroy:
            task_size = 30
        elif loner_type == "U-M" and level == 1 and self._destroy:
            task_size = 15
        else:
            if level % 2 == 0:
                task_size = 4
            else:
                task_size = 5
        return task_size

    def __handler_exit__(self,**kwargs):
        if self._constants.verbose:
            if self._total_loners - self._rescued["total"] == 0:
                self.__standard_output__("\n[Status] All reads succesfully paired") 
            else:
                self.__standard_output__("\n[Status] Couldn't find pairs for %d reads" %\
                    (self._total_loners - self._rescued["total"],))
        
        for qu in self._pause_qus:
            qu.close()

class ChaserClass(Process):
    def __init__(self,object parent_bam,object constants):
        super(ChaserClass,self).__init__()
        self._constants = constants

        head,root = os.path.split(parent_bam.filename)
        name,ext = os.path.splitext(root)
        self.basename = name

    def __sort_and_index__(self,path):
        pass

    def __get_loner_temp_path__(self,loner_type,level=0,unique=""):
        return "%s/%s_%s_lvl%d_typ%s%s.bam" %\
            (self._constants.temp_dir,self.basename,self.pid,level,loner_type,unique)

class PrimaryTask(ChaserClass):

    def __init__(self,object unsorted_paths,object outqu,object constants, int max_reads, dict child_pack):
        super(PrimaryTask,self).__init__(child_pack["parent_bam"],constants)
        
        self._unsorted_paths = unsorted_paths
        self._outqu = outqu
        self._child_pack = child_pack

        unsorted_object = pysam.AlignmentFile(unsorted_paths[0],"rb")
        self._references = unsorted_object.references
        self._header = unsorted_object.header
        unsorted_object.close()

        self._max_reads = max_reads

    def run(self):
        #speedup
        write_loner = self.__write_loner__
        #speedup

        loner_holder = self.__get_loner_holder__()
        loner_file_counts = Counter()

        loner_total = 0

        classifier_bins = self.__get_classifier_bins__()

        for read in self.__read_generator__():
            
            loner_type = self.__get_loner_type__(read,classifier_bins)
            loner_total += 1
            write_loner(loner_holder,loner_file_counts,loner_type,read)

        loner_paths = self.__get_paths_and_close__(loner_holder)

        for loner_type,paths in loner_paths.items():
            for path in paths:
                self.__start_matchmaker_task__([path],loner_type,-1) 

        del loner_paths,loner_file_counts,loner_holder

    def __start_matchmaker_task__(self,paths,loner_type,level):
        new_task = MatchMakerTask(paths,self._outqu,self._constants,
                                  loner_type,level,self._max_reads,
                                  self._child_pack)
        new_task.start()
        new_task.join()
        time.sleep(1)

    def __read_generator__(self):
        for path in self._unsorted_paths:
            self.__sort_and_index__(path)
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
            os.remove(path)
        gc.collect()

    def __get_classifier_bins__(self):
        crucial_bin_number = 30
        reference_len = len(self._references)

        if reference_len > crucial_bin_number:
            bins = []
            curcial_bins = list(np.digitize( range(crucial_bin_number) , range(0,crucial_bin_number,3)))

            remainder_bins = list(np.digitize(range(crucial_bin_number,reference_len),
                        range(crucial_bin_number,reference_len,50)) + 10)

            bins.extend(curcial_bins)
            bins.extend(remainder_bins)

            return bins
        else:
            return list(np.digitize( range(reference_len) , range(0,reference_len,2)))

    def __get_loner_type__(self,read,bins):
        if not read.is_unmapped and not read.mate_is_unmapped:
            if read.reference_id == read.next_reference_id:
                return "XX%d-%d" % (read.reference_id,read.next_reference_id,)
            else:
                class_bins = map(lambda x : bins[x],self.__order_reference_number__(read.reference_id,read.next_reference_id))
                return "%d-%d" % tuple(class_bins)
        elif read.is_unmapped and read.mate_is_unmapped:
            return "U-U"
        else:
            return "U-M"

    def __order_reference_number__(self,read_ref,mate_ref):
        if min((read_ref,mate_ref)) == read_ref:
            return read_ref,mate_ref
        else:
            return mate_ref,read_ref

    def __get_loner_holder__(self):
        holder = {}
        return holder

    def __get_loner_object__(self,loner_type,unique):
        path = self.__get_loner_temp_path__(loner_type,0,unique="-%d" % (unique,) )
        return pysam.AlignmentFile(path,"wb",header=self._header)

    def __write_loner__(self,loner_objects,loner_file_counts,loner_type,read):
        try:
            loner_path = loner_objects[loner_type][-1].filename
            if loner_file_counts[loner_path] == (self._max_reads-1): #minus one so we are always one less than max
                new_bam_object = self.__get_loner_object__(loner_type,len(loner_objects[loner_type]))
                loner_objects[loner_type].append(new_bam_object)
        except KeyError:
            new_bam_object = self.__get_loner_object__(loner_type,0)
            loner_objects[loner_type] = [new_bam_object]
        
        loner_objects[loner_type][-1].write(read)

    def __get_paths_and_close__(self,loner_holder):
        paths = {}
        for loner_type,loner_objects in loner_holder.items():
            paths[loner_type] = []
            for loner_object in loner_objects:
                paths[loner_type].append(loner_object.filename)
                loner_object.close()
        return paths

class MatchMakerTask(ChaserClass):

    def __init__(self,object loner_paths,object outqu,object constants,
        str loner_type,int level,int max_reads,object child_pack):

        super(MatchMakerTask,self).__init__(child_pack["parent_bam"],constants)
        
        self._max_reads = max_reads
        self._child_pack = child_pack
        self._loner_paths = loner_paths
        self._outqu = outqu
        self._loner_type = loner_type
        self._level = level
        self._header = self.__get_header__(self._loner_paths[0])
        self._leftover_paths = []

    def run(self):

        loner_path = self.__get_loner_temp_path__(self._loner_type,self._level+1)
        loner_object = self.__get_loner_object__(loner_path)

        cdef dict pairs = self.__first_pass__()
        cdef dict read_pairs = {}
        cdef dict send_pairs = {}
        cdef int rescued = 0
        cdef int written = 0

        for read in self.__read_generator__(second_pass=True):
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

                #we are handling pairs so we divide chunk by two.
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
            return_paths.append((self._level+1,loner_path))

        if len(self._leftover_paths) > 0:
            for leftover_path in self._leftover_paths:
                return_paths.append((self._level,leftover_path))

        loner_pack = MatchMakerPackage(loner_type=self._loner_type,results=return_paths,level=self._level+1,
                                       chaser_type="match_maker",rescued=rescued)
        self._outqu.put(loner_pack)

    def __get_header__(self,path):
        bam_object = pysam.AlignmentFile(path,"rb")
        header = bam_object.header
        bam_object.close()
        return header

    def __first_pass__(self):
        cdef dict pairs = {}
        for read in self.__read_generator__(second_pass=False):
            try:
                status = pairs[read.qname]
                pairs[read.qname] = True
            except KeyError:
                pairs[read.qname] = False
        return self.__prepare_for_second_pass__(pairs)

    def __prepare_for_second_pass__(self,pairs):
        return dict(filter( lambda tup : tup[1], iter(pairs.items()) ))

    def __read_generator__(self,second_pass=False):
        count = self._max_reads // len(self._loner_paths)
        count_limit = False
        for path in self._loner_paths:
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
                    leftover_object = self.__get_leftover_object__(len(self._leftover_paths))
                    for read in bam_iterator:
                        leftover_count += 1
                        leftover_object.write(read)
                    leftover_object.close()
                    
                    if leftover_count > 0:
                        self._leftover_paths.append(leftover_object.filename)
                    else:
                        os.remove(leftover_object.filename)
                os.remove(path)
            loner_object.close()
        gc.collect()

    def __launch_child_task__(self,pairs):
        child_pack = self._child_pack
        results = child_pack["chaser_task"].handle_pair_dict(pairs,self.pid)
        child_pack["queue"].put(parabam.core.Package(results=results))
        return len(pairs)

    def __get_leftover_object__(self,unique):
        path = self.__get_loner_temp_path__(self._loner_type,self._level,"-%d" % (unique,))
        leftover = pysam.AlignmentFile(path,"wb",header=self._header)
        return leftover

    def __get_loner_object__(self,path):
        loner_object = pysam.AlignmentFile(path,"wb",header=self._header)
        return loner_object

#.....happily ever after
