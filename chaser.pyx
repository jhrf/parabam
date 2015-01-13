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

from itertools import izip
from collections import Counter
from multiprocessing import Queue,Process
from parabam.core import Package

class ChaserPackage(Package):
    def __init__(self,name,results,destroy,level,source,chaser_type,total,sub_classifier):
        super(ChaserPackage,self).__init__(name,results,destroy)
        self.level = level
        self.chaser_type=chaser_type
        self.source = source
        self.total = total        
        self.sub_classifier = sub_classifier

class HandlerChaser(parabam.core.Handler):

    def __init__(self,object inqu,object mainqu,object const,object destroy_limit, object TaskClass):
        super(HandlerChaser,self).__init__(inqu,const,destroy_limit=destroy_limit,report=False)

        self._loner_types = ["both_mapped","both_unmapped","one_mapped"]
        self._sources = const.sources
        self._subset_types = const.subset_types
        
        self._mainqu = mainqu
        self._TaskClass = TaskClass
        
        self._loner_pyramid = self.__instalise_loner_pyramid__()

        self._periodic_interval = 5
        self._total_loners = 0
        self._rescued = Counter()

        self._match_tasks = []
        self._primary_tasks = []

    #borrowed from processor
    def __update_tasks__(self,active_tasks):
        terminated_procs = []
        for pr in active_tasks:
            if not pr.is_alive():
                pr.terminate()
                terminated_procs.append(pr)
        
        for pr in terminated_procs:
            active_tasks.remove(pr)

        if len(terminated_procs) > 0:
            #Only invoke the GC if there is the possibility of
            #collecting any memory.
            del terminated_procs
            gc.collect()

        return(len(active_tasks))

    def __instalise_loner_pyramid__(self):
        loner_pyramid = {}
        for source in self.const.sources:
            loner_pyramid[source] = {}
            for loner_type in self._loner_types:
                loner_pyramid[source][loner_type] = {}

        return loner_pyramid

    def __new_package_action__(self,new_package,**kwargs):
        if hasattr(new_package,"subset_type"):
            #This is a merge package
            self.__start_primary_task__(new_package)
        else:
            chaser_type = new_package.chaser_type
            if chaser_type == "primary":
                self.__handle_primary_task__(new_package)
            elif chaser_type == "match_maker":
                self.__handle_match_maker__(new_package)
            else:
                print "Error!"
 
    def __handle_match_maker__(self,new_package):
        source = new_package.source
        name = new_package.name
        level = new_package.level

        self._rescued["total"] += (new_package.total*2) 
        self._rescued[name] += (new_package.total*2) #this number is by pairs, to get read num we times 2

        pyramid = self._loner_pyramid[source][name]
        if level == -1:
            print "Handle The Terminal Case! The string collection is larger than the max size (5million?)"
        else:
            if new_package.results: #This ensures that there are leftover loners

                sub_pyramid = pyramid[new_package.sub_classifier]
                if level >= len(sub_pyramid):
                    sub_pyramid.append([])
                sub_pyramid[level].append(new_package.results)
                
    def __handle_primary_task__(self,new_package):
        source = new_package.source
        for loner_type,sub_classifier_dict in new_package.results.items():
            for sub_classifier,path in sub_classifier_dict.items():
                try: 
                    pyramid = self._loner_pyramid[source][loner_type][sub_classifier][new_package.level]
                except KeyError:
                    self._loner_pyramid[source][loner_type][sub_classifier] = [[]]
                    pyramid = self._loner_pyramid[source][loner_type][sub_classifier][new_package.level]

                pyramid.append(path)

    def __start_primary_task__(self,new_package):
        tasks_running = self.__update_tasks__(self._primary_tasks)
        if tasks_running <= 5:
            for loner_count,path in new_package.results:
                self._total_loners += loner_count
                new_task = PrimaryTask(path,self._inqu,self.const,new_package.source,self._loner_types)
                new_task.start()
                self._primary_tasks.append(new_task)
        else:
            self._inqu.put(new_package)

    def __periodic_action__(self,iterations):
        if iterations % 10 == 0 :
            sys.stdout.write(self.__format_out_str__())                
            sys.stdout.flush()

        for source,loner_dict in self._loner_pyramid.items():
            for loner_type,pyramid_dict in loner_dict.items():
                tasks_running = self.__update_tasks__(self._match_tasks)
                if tasks_running < 10:
                    for sub_classifier,pyramid in pyramid_dict.items():    
                        for i,sub_pyramid in enumerate(reversed(pyramid)):
                            level = len(pyramid) - (i+1)
                            task_size = 5
                            if len(sub_pyramid) >= task_size:
                                paths = [sub_pyramid.pop() for x in xrange(task_size)]
                                self.__start_match_maker__(paths,source,loner_type,level,sub_classifier)
                                break
        gc.collect()
        
    def __format_out_str__(self):
        induvidual_str = ""
        for i,loner_type in enumerate(self._loner_types):
            induvidual_str += "%i:%d " % (i,self._rescued[loner_type])

        return "\r %s Total: %d/%d Ratio: %.5f M-Task:%d P-Task:%d " %\
            (induvidual_str,self._rescued["total"],self._total_loners,
            float(self._rescued["total"]+1)/(self._total_loners+1),
            len(self._match_tasks),len(self._primary_tasks))

    def __get_task_size__(self,level,first_task,exponant):
        return (first_task // ( (level+1) ** exponant)) + 2

    def __start_match_maker__(self,paths,source,loner_type,level,sub_classifier):
        new_task = MatchMakerTask(paths,self._inqu,self.const,source,loner_type,level,sub_classifier,
                    {"queue":self._mainqu,"const":self.const,
                     "task_args":[source],"TaskClass":self._TaskClass})
        new_task.start()
        self._match_tasks.append(new_task)

    def __handler_exit__(self,**kwargs):
        print "Exiting"

# class EndOfLonerProcessing(object):

#     def __init__(self,pyramid,object const,TaskClass,task_args,queue):

#     def run(self):
#         self.__wait_for_all_tasks__()
#         self.__condense_existing_pyramid__(self._pyramid)

#     def __condense_existing_pyramid__(self):

#         for source,loner_dict in self._pyramid.items():
#             for loner_type,pyramid in loner_dict.items():
#                 if len(self._match_tasks) < 5:
#                     for i,sub_pyramid in enumerate(pyramid):
#                         paths = [sub_pyramid.pop() for x in xrange(task_size)]
#                         new_task = MatchMakerTask(paths,self._inqu,self.const,source,loner_type,level,
#                                     {"queue":self._mainqu,"const":self.const,
#                                      "task_args":[source],"TaskClass":self._TaskClass})
#                         new_task.start()
#                         new_task.join()

#                         try:
#                             results = self.queue.get(100)
#                             #put this into the pyramid
#                         except Queue2.Empty:
#                             pass

#     def __top_level_comparison__(self):

#     def __finish_subset_analysis__(self):

#     def __merge_chaser_to_existing__(self):

class ChaserClass(Process):
    def __init__(self,object const,str source):
        super(ChaserClass,self).__init__()
        self.const = const
        self._source = source

    def __sort_and_index__(self,path):
        temp_path = "%s/sort_temp_%d_%s" % (self.const.temp_dir,time.time(),self.pid)
        pysam.sort(path,temp_path)
        temp_path += ".bam"

        os.remove(path)
        os.rename(temp_path,path)

        pysam.index(path)

    def __get_loner_temp_path__(self,loner_type,sub_classifier,level=0):
        return "%s/%s_%s_%s_%d_level_%d_%s.bam" %\
            (self.const.temp_dir,self._source,loner_type,self.pid,time.time(),level,sub_classifier)

class PrimaryTask(ChaserClass):

    def __init__(self,object unsorted_path,object outqu,object const,str source,list loner_types):
        super(PrimaryTask,self).__init__(const,source)
        
        self._unsorted_path = unsorted_path
        self._outqu = outqu
        self._loner_types = loner_types

    def run(self):
        #speedup
        write_loner = self.__write_loner__
        #speedup

        unsorted_path = self._unsorted_path

        self.__sort_and_index__(unsorted_path)
        unsorted_object = pysam.AlignmentFile(unsorted_path,"rb")

        loner_holder = self.__get_loner_holder__()
        loner_counts = Counter()
        loner_total = 0

        classifier_bins = np.digitize(range(len(unsorted_object.references)),range(0,len(unsorted_object.references),10))

        for read in self.__read_generator__(unsorted_object):
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

        loner_paths = self.__get_paths_and_close__(loner_holder)

        unsorted_object.close()
        os.remove(unsorted_path)
        os.remove(unsorted_path+".bai")

        loner_pack = ChaserPackage(name="all",results=loner_paths,destroy=False,level=0,source=self._source,
            chaser_type="primary",total=loner_total,sub_classifier=None)
        self._outqu.put(loner_pack)

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
            return "u-%d" % (bins[self.__order_reference_number__(read.reference_id,read.next_reference_id)[1]],)
        else:
            return "e-e"

    def __order_reference_number__(self,read_ref,mate_ref):
        if min((read_ref,mate_ref)) == read_ref:
            return read_ref,mate_ref
        else:
            return mate_ref,read_ref

    def __read_generator__(self,alig_file_obj):
        for read in alig_file_obj.fetch(until_eof=True):
            yield read

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
        str loner_type,int level,str sub_classifier,object child_package):

        super(MatchMakerTask,self).__init__(const,source)
        
        self._child_package = child_package
        self._loner_paths = loner_paths
        self._outqu = outqu
        self._loner_type = loner_type
        self._level = level
        self._sub_classifier = sub_classifier

    def run(self):

        from pprint import pprint as ppr 

        loner_path = self.__get_loner_temp_path__(self._loner_type,self._sub_classifier,self._level+1)
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
                if len(send_pairs) >= (self.const.chunk//2):
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

        if written == 0:
            os.remove(loner_path)
            loner_path = None

        loner_pack = ChaserPackage(name=self._loner_type,results=loner_path,destroy=False,level=self._level+1,
                source=self._source,chaser_type="match_maker",total=rescued,sub_classifier=self._sub_classifier)
        self._outqu.put(loner_pack)

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
        template = pysam.AlignmentFile(self._loner_paths[0],"rb")
        loner_object = pysam.AlignmentFile(path,"wb",template=template)
        template.close()
        return loner_object

#.....happily ever after
