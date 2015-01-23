#Once upon a time
import pysam
import pdb
import parabam
import time
import sys
import os
import gc

import Queue as Queue2
import numpy as np

from random import randint
from itertools import izip
from collections import Counter
from multiprocessing import Queue,Process
from parabam.core import Package,CorePackage

class ChaserPackage(Package):
    def __init__(self,name,results,destroy,level,source,chaser_type,total,sub_classifier,order,processing=True):
        super(ChaserPackage,self).__init__(name,results,destroy)
        self.level = level
        self.chaser_type=chaser_type
        self.source = source
        self.total = total        
        self.sub_classifier = sub_classifier
        self.order = order
        self.processing = processing

class HandlerChaser(parabam.core.Handler):

    def __init__(self,object inqu,object mainqu,object pause_qu,object const,object destroy_limit, object TaskClass):
        super(HandlerChaser,self).__init__(inqu,const,destroy_limit=destroy_limit,report=False)

        self._loner_types = ["both_mapped","both_unmapped","one_mapped"]
        self._sources = const.sources
        self._subset_types = const.subset_types
        
        self._pause_qu = pause_qu
        self._mainqu = mainqu
        self._TaskClass = TaskClass
        
        self._loner_pyramid = self.__instalise_loner_pyramid__()
        self._pyramid_idle_counts = Counter()

        self._periodic_interval = 5
        self._total_loners = 0
        self._rescued = Counter()

        self._tasks = []
        self._primary_store = self.__instalise_primary_store__()

        self._processing = True

    #START -- BORROWED FROM PROCESOR
    #Need to refactor here to stop code duplication
    def __wait_for_tasks__(self,list active_tasks,int max_tasks=8):
        update_tasks = self.__update_tasks__ #optimising alias
        update_tasks(active_tasks)
        cdef int currently_active = len(active_tasks)
        self._active_count = currently_active

        if max_tasks > currently_active:
            return

        self._pause_qu.put(True)
        while(max_tasks <= currently_active and currently_active > 0):
            update_tasks(active_tasks)
            currently_active = len(active_tasks)
            time.sleep(1)
        self._pause_qu.put(False)

        return 

    def __instalise_primary_store__(self):
        paths = {}
        for source in self._sources:
            paths[source] = []
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

    #END -- BORROWED FROM PROCESOR

    def __instalise_loner_pyramid__(self):
        loner_pyramid = {}
        for source in self.const.sources:
            loner_pyramid[source] = {}
            for loner_type in self._loner_types:
                loner_pyramid[source][loner_type] = {}

        return loner_pyramid

    def __new_package_action__(self,new_package,**kwargs):

        chaser_type = new_package.chaser_type
        if chaser_type == "origin":
            self.__handle_origin_task__(new_package)
        elif chaser_type == "primary":
            self.__handle_primary_task__(new_package)
        elif chaser_type == "match_maker":
            self.__handle_match_maker__(new_package)
        else:
            print "Error!"

    def __handle_primary_task__(self,new_package):
        source = new_package.source
        for loner_type,sub_classifier_dict in new_package.results.items():
            for sub_classifier,path in sub_classifier_dict.items():
                try: 
                    pyramid = self._loner_pyramid[source][loner_type][sub_classifier][new_package.level]
                except KeyError:
                    self._loner_pyramid[source][loner_type][sub_classifier] = [[]]
                    pyramid = self._loner_pyramid[source][loner_type][sub_classifier][new_package.level]
                pyramid.append((new_package.order,path))
                self._pyramid_idle_counts[(loner_type,sub_classifier)] = 0

    def __handle_origin_task__(self,new_package):
        self._processing = new_package.processing

        for loner_count,path in new_package.results:
            if loner_count == 0:
                os.remove(path)
            else:
                self._total_loners += loner_count
                self._primary_store[new_package.source].append((new_package.order,path))

        self.__test_primary_tasks__()

    def __handle_match_maker__(self,new_package):
        source = new_package.source
        name = new_package.name
        level = new_package.level

        self._rescued["total"] += (new_package.total*2) 
        self._rescued[name] += (new_package.total*2) #this number is by pairs, to get read num we times 2

        if new_package.results: #This ensures that there are leftover loners
            pyramid = self._loner_pyramid[source][name][new_package.sub_classifier]
            if level >= len(pyramid):
                pyramid.append([])
            if len(new_package.results) == 1: 
                pyramid[level].append((new_package.order,new_package.results[0]))
            else:
                for loner_path in new_package.results:
                    random_insert = randint(0,len(pyramid[level]))
                    pyramid[level].insert(random_insert,(new_package.order,loner_path))

    def __test_primary_tasks__(self,required_paths=10):
        task_sent = []
        for source,paths in self._primary_store.items():
            if len(paths) >= required_paths:
                self.__start_primary_task__(source)
                task_sent.append(source)

        for source in task_sent:
            self._primary_store[source] = []
        gc.collect()
          
    def __start_primary_task__(self,source):        
        self.__wait_for_tasks__(self._tasks) #wait for existing tasks

        orders,paths = zip(*self._primary_store[source])

        new_task = PrimaryTask(paths,self._inqu,self.const,source,self._loner_types,max(orders))
        new_task.start()
        self._tasks.append(new_task)

    def __start_matchmaker_task__(self,paths,source,loner_type,level,sub_classifier,max_order):
        self.__wait_for_tasks__(self._tasks)

        new_task = MatchMakerTask(paths,self._inqu,self.const,source,
                    loner_type,level,sub_classifier,max_order,
                    {"queue":self._mainqu,"const":self.const,
                     "task_args":[source],"TaskClass":self._TaskClass})
        new_task.start()

        self._tasks.append(new_task)
        return max_order

    def __periodic_action__(self,iterations):
        if self._processing:
            idle_threshold = 300
            self.__test_primary_tasks__()
        else:
            idle_threshold = 5
            self.__test_primary_tasks__(required_paths=1)

        if iterations % 10 == 0 :
            sys.stdout.write(self.__format_out_str__())                
            sys.stdout.flush()

        pyramid_idle_counts = self._pyramid_idle_counts 
        empty = True
        for source,loner_dict in self._loner_pyramid.items():
            for loner_type,pyramid_dict in loner_dict.items():
                    for sub_classifier,pyramid in pyramid_dict.items():
                        for i,sub_pyramid in enumerate(reversed(pyramid)):
                            level = len(pyramid) - (i+1)
                            task_size = self.__get_task_size__(level,sub_classifier)
                            has_mix,max_order = self.__order_mix__(task_size,sub_pyramid)
                            if len(sub_pyramid) >= task_size and has_mix:
                                empty = False
                                paths = [sub_pyramid.pop()[1] for x in xrange(task_size)]
                                self.__start_matchmaker_task__(paths,source,loner_type,level,sub_classifier,max_order)
                                pyramid_idle_counts[(loner_type,sub_classifier)] = 0
                                break

                            elif len(sub_pyramid) > 0:
                                empty = False
                                pyramid_idle_counts[(loner_type,sub_classifier)] += 1
                                if len(sub_pyramid) > task_size:
                                    self.__cycle_paths__(sub_pyramid)

                        if pyramid_idle_counts[(loner_type,sub_classifier)] >= idle_threshold:
                            idle_success,idle_paths = self.__idle_routine__(pyramid,sub_classifier,loner_type,source)
                            pyramid_idle_counts[(loner_type,sub_classifier)] = 0

                            #delete after idle success or when processing has finished and idle fails
                            if self.__clear_subpyramid__(idle_success,self._processing):  
                                for delete_tup in idle_paths:
                                    for sub_pyramid in pyramid:
                                        if delete_tup in sub_pyramid:
                                            sub_pyramid.remove(delete_tup)
                                            break
        
        self.__update_tasks__(self._tasks)              
        if empty and not self._processing and len(self._tasks) == 0:
            print "\n[Status] Couldn't find pairs for %d reads" % (self._total_loners - self._rescued["total"],)
            self.__send_final_kill_signal__()

        gc.collect()

    def __cycle_paths__(self,sub_pyramid):
        path = sub_pyramid.pop()
        sub_pyramid.append(path)

    def __order_mix__(self,task_size,sub_pyramid):
        seen = []
        for i,(order,path) in enumerate(sub_pyramid):
            if i == task_size:
                break
            seen.append(order)
        seen_set = set(seen)
        if len(seen_set) > 1:
            return True,max(seen_set)
        else:
            return False,-1

    def __clear_subpyramid__(self,success,processing):
        if success:
            #If the idle task was sent, clear pyramid
            return True
        elif not processing and not success:
            running = self.__update_tasks__(self._tasks)
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

    def __send_final_kill_signal__(self):
        kill_package = CorePackage(name="",results={},destroy=True,curproc=0)
        self._mainqu.put(kill_package)

    def __idle_routine__(self,pyramid,sub_classifier,loner_type,source):
        paths = []
        orders = []
        task_size = 3
        max_level = 0

        for i,sub_pyramid in enumerate(pyramid):
            max_level = i
            for order,path in sub_pyramid:
                if len(orders) == 0 or any(map(lambda comp: not comp == order,orders)):
                    paths.append(path)
                    orders.append(order)
                    if len(paths) == task_size:
                        break
        
        if (len(paths) > 1) or (len(paths) > 0 and max_level == 0): #level 0 may not be unique
            order = self.__start_matchmaker_task__(paths,source,loner_type,max_level,sub_classifier,max(orders))
            return True,zip(orders,paths)

        return False,zip(orders,paths)
       
    def __get_task_size__(self,level,sub_classifier):
        if sub_classifier == "0-0":
            return 2
        elif sub_classifier == "u-m" and level == 0:
            return 20
        else:
            return 6

    def __format_out_str__(self):
        induvidual_str = ""
        for i,loner_type in enumerate(self._loner_types):
            induvidual_str += "%i:%d " % (i,self._rescued[loner_type])

        return "\r %s Total: %d/%d Ratio: %.5f Tasks:%d " %\
            (induvidual_str,self._rescued["total"],self._total_loners,
            float(self._rescued["total"]+1)/(self._total_loners+1),
            len(self._tasks))

    def __handler_exit__(self,**kwargs):
        pass

class ChaserClass(Process):
    def __init__(self,object const,str source):
        super(ChaserClass,self).__init__()
        self.const = const
        self._source = source

    def __sort_and_index__(self,path):
        pass

    def __get_loner_temp_path__(self,loner_type,sub_classifier,level=0,unique=""):
        return "%s/%s_%s_%s_%d_level_%d_%s%s.bam" %\
            (self.const.temp_dir,self._source,loner_type,self.pid,time.time(),level,sub_classifier,unique)

class PrimaryTask(ChaserClass):

    def __init__(self,object unsorted_paths,object outqu,object const,str source,list loner_types,int order):
        super(PrimaryTask,self).__init__(const,source)
        
        self._order = order
        self._unsorted_paths = unsorted_paths
        self._unsorted_path = unsorted_paths[0]
        self._outqu = outqu
        self._loner_types = loner_types

    def run(self):
        #speedup
        write_loner = self.__write_loner__
        #speedup

        unsorted_path = self._unsorted_path
        unsorted_object = pysam.AlignmentFile(unsorted_path,"rb")

        loner_holder = self.__get_loner_holder__()
        loner_counts = Counter()
        loner_total = 0

        classifier_bins = self.__get_classifier_bins__(unsorted_object)

        for read in self.__read_generator__():
            if read.is_unmapped and read.mate_is_unmapped:
                current_type = "both_unmapped"
            elif not read.is_unmapped and not read.mate_is_unmapped:
                current_type = "both_mapped"
            else:
                current_type = "one_mapped"

            sub_classifier = self.__get_sub_classifier__(current_type,read,classifier_bins)
            write_loner(loner_holder,current_type,sub_classifier,unsorted_object,read)

            loner_counts[current_type] += 1
            loner_total += 1

        unsorted_object.close()
        loner_paths = self.__get_paths_and_close__(loner_holder)

        loner_pack = ChaserPackage(name="all",results=loner_paths,destroy=False,level=0,source=self._source,
            chaser_type="primary",total=loner_total,sub_classifier=None,order=self._order)
        self._outqu.put(loner_pack)

    def __read_generator__(self,second_pass=False):
        for path in self._loner_paths:
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
            if second_pass:
                os.remove(path)
        gc.collect()

    def __read_generator__(self):
        for path in self._unsorted_paths:
            self.__sort_and_index__(path)
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
            os.remove(path)
            #os.remove(path+".bai")
        gc.collect()

    def __get_classifier_bins__(self,unsorted_object):
        crucial_bin_number = 30
        if len(unsorted_object.references) > crucial_bin_number:
            bins = []
            curcial_bins = list(np.digitize( range(crucial_bin_number) , range(0,crucial_bin_number,3)))

            remainder_bins = list(np.digitize(range(crucial_bin_number,len(unsorted_object.references)),
                        range(crucial_bin_number,len(unsorted_object.references),50)) + 10)

            bins.extend(curcial_bins)
            bins.extend(remainder_bins)

            return bins
        else:
            return list(np.digitize( range(len(unsorted_object.references)) , range(0,len(unsorted_object.references),2)))

    def __get_sub_classifier__(self,loner_type,read,bins):
        if loner_type == "both_mapped":
            if read.reference_id == read.next_reference_id:
                return "0-0"
            else:
                class_bins = map(lambda x : bins[x],self.__order_reference_number__(read.reference_id,read.next_reference_id))
                return "%d-%d" % tuple(class_bins)
        elif loner_type == "both_unmapped":
            return "u-u"
        elif loner_type == "one_mapped":
            return "u-m" #mate reference not loaded properly
        else:
            return "e-e"

    def __order_reference_number__(self,read_ref,mate_ref):
        if min((read_ref,mate_ref)) == read_ref:
            return read_ref,mate_ref
        else:
            return mate_ref,read_ref

    def __get_loner_holder__(self):
        holder = {}
        for loner_type in self._loner_types:
            holder[loner_type] = {}
        return holder

    def __get_loner_object__(self,loner_type,sub_classifier,unsorted_object):
        path = self.__get_loner_temp_path__(loner_type,sub_classifier,0)
        return pysam.AlignmentFile(path,"wb",template=unsorted_object)

    def __write_loner__(self,loner_objects,loner_type,sub_classifier,unsorted_object,read):
        try:
            loner_objects[loner_type][sub_classifier].write(read)
        except KeyError:
            new_bam_object = self.__get_loner_object__(loner_type,sub_classifier,unsorted_object)
            loner_objects[loner_type][sub_classifier] = new_bam_object
            loner_objects[loner_type][sub_classifier].write(read)

    def __get_paths_and_close__(self,loner_holder):
        paths = {}
        for loner_type,sub_classifier_dict in loner_holder.items():
            paths[loner_type] = {}
            for sub_classifier,loner_object in sub_classifier_dict.items():
                paths[loner_type][sub_classifier] = loner_object.filename
                loner_object.close()
        return paths

class MatchMakerTask(ChaserClass):

    def __init__(self,object loner_paths,object outqu,object const,str source,
        str loner_type,int level,str sub_classifier,int max_order,object child_package):

        super(MatchMakerTask,self).__init__(const,source)
        
        self._child_package = child_package
        self._loner_paths = loner_paths
        self._outqu = outqu
        self._loner_type = loner_type
        self._level = level
        self._sub_classifier = sub_classifier
        self._max_order = max_order
        self._header = self.__get_header__(self._loner_paths[0])

    def run(self):
        loner_paths = [] 
        loner_objects = []
        self.__new_loner_setup__(loner_paths,loner_objects)

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
                if len(send_pairs) >= (self.const.chunk//2):
                    rescued += self.__launch_child_task__(send_pairs)
                    send_pairs = {}
                    gc.collect()

            except KeyError:
                written += 1
                if written % 2500000 == 0:
                    self.__new_loner_setup__(loner_paths,loner_objects)           
                loner_objects[-1].write(read)

        loner_objects[-1].close()

        if len(send_pairs) > 0:
            rescued += self.__launch_child_task__(send_pairs)

        del read_pairs   
        del pairs

        gc.collect()

        if written == 0:
            os.remove(loner_paths[-1])
            loner_path = None

        loner_pack = ChaserPackage(name=self._loner_type,results=loner_paths,destroy=False,level=self._level+1,
                source=self._source,chaser_type="match_maker",total=rescued,sub_classifier=self._sub_classifier,order=self._max_order)
        self._outqu.put(loner_pack)

    def __get_header__(self,path):
        bam_object = pysam.AlignmentFile(path,"rb")
        header = bam_object.header
        bam_object.close()
        return header

    def __new_loner_setup__(self,loner_paths,loner_objects):
        unique = ""
        if len(loner_objects) > 0:
            loner_objects[-1].close()
            unique = "_%d" % (len(loner_objects),)

        loner_paths.append(self.__get_loner_temp_path__(self._loner_type,self._sub_classifier,self._level+1,unique=unique))
        loner_objects.append(self.__get_loner_object__(loner_paths[-1]))

    def __remove_paths__(self):
        for path in self._loner_paths:
            os.remove(path)

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
        for path in self._loner_paths:
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
            if second_pass:
                os.remove(path)
        gc.collect()

    def __launch_child_task__(self,pairs):
        child_pack = self._child_package
        args = [pairs,child_pack["queue"],0,False,child_pack["const"]]
        args.extend(child_pack["task_args"])
        task = child_pack["TaskClass"](*args)
        task.start()
        task.join()
        return len(pairs)

    def __get_loner_object__(self,path):
        loner_object = pysam.AlignmentFile(path,"wb",header=self._header)
        return loner_object

#.....happily ever after
