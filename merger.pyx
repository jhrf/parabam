#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os,argparse
import pysam
import time
import parabam
import gzip

import time

import Queue as Queue2

from collections import Counter
from multiprocessing import Queue,Process
from itertools import izip

class MergePackage(parabam.core.Package):
    def __init__(self,object results,str subset_type,str source,object destroy):
        super(MergePackage,self).__init__(results,destroy)
        self.subset_type = subset_type
        self.source = source

class Handler(parabam.core.Handler):
    def __init__(self,object inqu, object const,int destroy_limit=1):
        super(Handler,self).__init__(inqu,const,report=False,destroy_limit=destroy_limit)
        
        self._sources = const.sources
        self._user_subsets = list(const.user_subsets)
        self._out_file_objects = self.__get_out_file_objects__()
        self._merged = 0

    def __get_out_file_objects__(self):
        file_objects = {}
        output_paths = self.const.output_paths
        for src in self._sources:
            master = pysam.Samfile(self.const.master_file_path[src],"rb")
            file_objects[src] = {}
            for subset in self._user_subsets:
                merge_type = self.__get_subset_merge_type__(output_paths[src][subset])
                cur_file_obj = self.__get_file_for_merge_type__(merge_type,output_paths[src][subset],master)
                file_objects[src][subset] = cur_file_obj
            master.close()
        return file_objects

    def __get_file_for_merge_type__(self,merge,path,master):
        if merge == "gzip":
            return gzip.open(path,"wb")
        elif merge == "txt":
            return open(path,"w")
        return pysam.AlignmentFile(path,"wb",template=master)

    def __get_subset_merge_type__(self,path):
        root,extension = os.path.splitext(path)
        if extension == ".gzip" or extension == ".gz":
            return "gzip"
        elif extension == ".txt":
            return "txt"
        return "pysam"

    def __periodic_action__(self,iterations):
        pass

    def __new_package_action__(self,new_package,**kwargs):
        #Handle the result. Result will always be of type parabam.support.MergePackage
        subset_type = new_package.subset_type
        source = new_package.source

        for merge_count,merge_path in new_package.results:
            try:
                if merge_count > 0:
                    for item in self.__get_entries_from_file__(merge_path,subset_type):
                        self._out_file_objects[source][subset_type].write(item)
                os.remove(merge_path)

            except IOError:
                print "FAILURE ON THIS PATH %s" % (merge_path,)

                print "\nTEMP CONTENTS"
                sys.stdout.write("\n".join(os.listdir(self.const.temp_dir)))
                print "--"
                print "Printing newpackage.results"
                print new_package.results
                print "--"

                raise
                
            self._merged += 1       

    def __get_entries_from_file__(self,path,subset):
        merge_type = self.__get_subset_merge_type__(path)
        if merge_type == "pysam":
            file_object = pysam.AlignmentFile(path,"rb")
            for read in file_object.fetch(until_eof=True):
                yield read
        elif merge_type == "gzip":
            file_object = gzip.open(path,"rb")
            for item in file_object:
                yield item
        elif merge_type == "txt":
            file_object = open(path,"r")
            for item in file_object:
                yield item
        else:
            print "[Error] Unrecognised merge type"
            return
        file_object.close()

    def __get_unique_temp_path__(self,temp_type,temp_dir="."):
        return "%s/%sTEMP%d.bam" % (temp_dir,temp_type,int(time.time()),)

    def __close_all_out_files__(self):
        for source,file_dict in self._out_file_objects.items():
            for subset,file_obj in file_dict.items():
                file_obj.close()
        del self._out_file_objects

    def __handler_exit__(self,**kwargs):
        self.__close_all_out_files__()
        
#...happily ever after
