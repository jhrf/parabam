#Once upon a time
import pysam
import pdb
import parabam
import time
import sys
import os
import gc

from itertools import izip
from collections import Counter
from multiprocessing import Queue,Process
from parabam.core import Package

class ChaserPackage(Package):
    def __init__(self,name,results,destroy,level,source,chaser_type,total):
        super(ChaserPackage,self).__init__(name,results,destroy)
        self.level = level
        self.chaser_type=chaser_type
        self.source = source
        self.total = total        

class HandlerChaser(parabam.core.Handler):

    def __init__(self,object inqu,object mainqu,object const,object destroy_limit, object TaskClass):
        super(HandlerChaser,self).__init__(inqu,const,destroy_limit=destroy_limit,report=False)

        self._loner_types = ["both_mapped","both_unmapped","one_mapped"]
        self._sources = const.sources
        self._subset_types = const.subset_types
        
        self._mainqu = mainqu
        self._TaskClass = TaskClass
        
        self._loner_pyramid = self.__instalise_loner_pyramid__()

        self._periodic_interval = 10
        self._total_loners = 0
        self._rescued = 0

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

    def __instalise_loner_pyramid__(self):
        loner_pyramid = {}
        for source in self.const.sources:
            loner_pyramid[source] = {}
            for loner_type in self._loner_types:
                loner_pyramid[source][loner_type] = [[]]
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

        self._rescued += (new_package.total*2) #this number is by pairs, to get read num we times 2

        pyramid = self._loner_pyramid[source][name]
        if level == -1:
            print "Handle The Terminal Case! The string collection is larger than the max size (5million?)"
        else:
            if level >= len(pyramid):
                pyramid.append([])
            pyramid[level].append(new_package.results)

    def __handle_primary_task__(self,new_package):
        source = new_package.source
        for loner_type,loner_path in new_package.results.items():
            self._loner_pyramid[source][loner_type][0].append(loner_path)

    def __start_primary_task__(self,new_package):
        self.__update_tasks__(self._primary_tasks)
        if len(self._primary_tasks) < 10:
            for loner_count,path in new_package.results:
                self._total_loners += loner_count
                new_task = PrimaryTask(path,self._inqu,self.const,new_package.source,self._loner_types)
                new_task.start()
                self._primary_tasks.append(new_task)
        else:
            self._inqu.put(new_package)

    def __periodic_action__(self,iterations):
        if iterations % 30 ==0 :
            sys.stdout.write("\rRescued: %d Total Loners: %d Ratio: %.5f Match_Tasks: %d Primary_Tasks: %d" % \
                (self._rescued,self._total_loners,float(self._rescued+1)/(self._total_loners+1),len(self._match_tasks),len(self._primary_tasks)))
            sys.stdout.flush()

        task_size = 8

        for source,loner_dict in self._loner_pyramid.items():
            for loner_type,pyramid in loner_dict.items():
                self.__update_tasks__(self._match_tasks)
                if len(self._match_tasks) < 10:
                    for i,sub_pyramid in enumerate(reversed(pyramid)):
                        level = len(pyramid) - (i+1)
                        if len(sub_pyramid) >= task_size:
                            paths = [sub_pyramid.pop() for x in  xrange(task_size)]
                            new_task = MatchMakerTask(paths,self._inqu,self.const,source,loner_type,level,
                                        {"queue":self._mainqu,"const":self.const,
                                         "task_args":[source],"TaskClass":self._TaskClass})
                            new_task.start()
                            self._match_tasks.append(new_task)
                            break
   
    def __handler_exit__(self,**kwargs):
        print "Exiting..."

class ChaserClass(Process):
    def __init__(self):
        super(ChaserClass,self).__init__()

    def __sort_and_index__(self,path):
        temp_path = "%s/sort_temp_%d_%s" % (self.const.temp_dir,time.time(),self.pid)
        pysam.sort(path,temp_path)
        temp_path += ".bam"

        os.remove(path)
        os.rename(temp_path,path)

        pysam.index(path)

    def __get_loner_temp__(self,loner_type,level=0):
        return "%s/%s_%s_%d_level%d.bam" % (self.const.temp_dir,loner_type,self.pid,time.time(),level,)

class PrimaryTask(ChaserClass):

    def __init__(self,object unsorted_path,object outqu,object const,str source,list loner_types):
        super(PrimaryTask,self).__init__()
        
        self._unsorted_path = unsorted_path
        self._outqu = outqu
        self._source = source
        self.const = const
        self._loner_types = loner_types

    def run(self):
        #speedup
        write_loner = self.__write_loner_type__
        #speedup

        unsorted_path = self._unsorted_path

        self.__sort_and_index__(unsorted_path)
        unsorted_object = pysam.AlignmentFile(unsorted_path,"rb")

        loner_types = self._loner_types
        loner_paths = self.__create_loner_paths__(loner_types)
        loner_objects = self.__create_loner_objects__(loner_paths,unsorted_object)
        loner_counts = Counter()
        loner_total = 0

        for read in self.__read_generator__(unsorted_object):
            current_type = ""
            if read.is_unmapped and read.mate_is_unmapped:
                current_type = "both_unmapped"
            elif not read.is_unmapped and not read.mate_is_unmapped:
                current_type = "both_mapped"
            else:
                current_type = "one_mapped"

            loner_counts[current_type] += 1
            loner_total += 1
            write_loner(loner_objects,current_type,read)

        self.__close_loner_objects__(loner_objects)

        for loner_type in loner_types:
            #remove empties
            if loner_counts[loner_type] == 0:
                os.remove(loner_paths[loner_type])
                del loner_paths[loner_type]

        os.remove(unsorted_path)
        os.remove(unsorted_path+".bai")

        loner_pack = ChaserPackage(name="all",results=loner_paths,destroy=False,level=0,
            source=self._source,chaser_type="primary",total=loner_total)
        self._outqu.put(loner_pack)

    def __read_generator__(self,alig_file_obj):
        for read in alig_file_obj.fetch(until_eof=True):
            yield read

    def __create_loner_objects__(self,loner_paths,unsorted_object):
        objects = {}
        for loner_type,path in loner_paths.items():
            objects[loner_type] = pysam.AlignmentFile(path,"wb",template=unsorted_object)
        return objects

    def __create_loner_paths__(self,loner_types):
        paths = {}
        for loner_type in loner_types:
            paths[loner_type] = self.__get_loner_temp__(loner_type,0)
        return paths

    def __write_loner_type__(self,loner_objects,loner_type,read):
        loner_objects[loner_type].write(read)

    def __close_loner_objects__(self,loner_objects):
        for loner_type,loner_object in loner_objects.items():
            loner_object.close()

class MatchMakerTask(ChaserClass):

    def __init__(self,object loner_paths,object outqu,object const,str source,str loner_type,int level,object child_package):
        super(MatchMakerTask,self).__init__()
        
        self._child_package = child_package
        self._loner_paths = loner_paths
        self._outqu = outqu
        self._source = source
        self._loner_type = loner_type
        self._level = level
        self.const = const

    def run(self):

        from pprint import pprint as ppr 

        loner_path = self.__get_loner_temp__(self._loner_type,self._level+1)
        loner_object = self.__get_loner_object__(loner_path)

        cdef dict pairs = self.__first_pass__()
        cdef dict read_pairs = {}
        cdef dict send_pairs = {}
        cdef int rescued = 0

        for read in self.__read_generator__():
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
                loner_object.write(read)

        loner_object.close()

        if len(send_pairs) > 0:
            rescued += self.__launch_child_task__(send_pairs)

        del read_pairs   
        del pairs

        self.__remove_paths__()
        gc.collect()

        #self.__sort_and_index__(loner_path)
        loner_pack = ChaserPackage(name=self._loner_type,results=loner_path,destroy=False,level=self._level+1,
            source=self._source,chaser_type="match_maker",total=rescued)
        self._outqu.put(loner_pack)

    def __remove_paths__(self):
        for path in self._loner_paths:
            os.remove(path)
            #os.remove(path+".bai")

    def __first_pass__(self):
        cdef dict pairs = {}
        for read in self.__read_generator__():
            try:
                status = pairs[read.qname]
                pairs[read.qname] = True
            except KeyError:
                pairs[read.qname] = False
        return self.__prepare_for_second_pass__(pairs)

    def __prepare_for_second_pass__(self,pairs):
        return dict(filter( lambda tup : tup[1], iter(pairs.items()) ))

    def __read_generator__(self):
        for path in self._loner_paths:
            loner_object = pysam.AlignmentFile(path,"rb")
            for read in loner_object.fetch(until_eof=True):
                yield read
            loner_object.close()
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