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
    def __init__(self,object results,str subset_type):
        super(MergePackage,self).__init__(results)
        self.subset_type = subset_type

class Handler(parabam.core.Handler):

    def __init__(self,object parent_bam, object output_paths,object inqu,
                      object const,object pause_qu,dict out_qu_dict):

        super(Handler,self).__init__(parent_bam,output_paths,inqu,const,pause_qu,
                                        out_qu_dict,report=False)
        
        self._user_subsets = list(const.user_subsets)
        self._out_file_objects = self.__get_out_file_objects__()
        self._merged = 0

    def __get_out_file_objects__(self):
        file_objects = {}
        output_paths = self._output_paths

        file_objects = {}
        for subset in self._user_subsets:
            output_path = output_paths[self._parent_bam.filename][subset]
            merge_type = self.__get_subset_merge_type__(output_path)
            cur_file_obj = self.__get_file_for_merge_type__(merge_type,output_path)
            file_objects[subset] = cur_file_obj

        return file_objects

    def __get_file_for_merge_type__(self,merge,path):
        if merge == "gzip":
            return gzip.open(path,"wb")
        elif merge == "txt":
            return open(path,"w")
        return pysam.AlignmentFile(path,"wb",header=self._parent_bam.header)

    def __get_subset_merge_type__(self,path):
        root,extension = os.path.splitext(path)
        if extension == ".gzip" or extension == ".gz":
            return "gzip"
        elif extension == ".txt":
            return "txt"
        return "pysam"

    def __periodic_action__(self,iterations):
        if self._destroy:
            try:
                pack = self._inqu.get(False,20)
                self._inqu.put(pack)
            except Queue2.Empty:
                self._finished = True

    def __new_package_action__(self,new_package,**kwargs):
        #Handle the result. Result will always be of type parabam.support.MergePackage
        subset_type = new_package.subset_type
        self._pause_qu.put(True)
        for merge_count,merge_path in new_package.results:
            try:
                if merge_count > 0:
                    for item in self.__get_entries_from_file__(merge_path,subset_type):
                        self._out_file_objects[subset_type].write(item)
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
        time.sleep(1)
        self._pause_qu.put(False)

    def __get_entries_from_file__(self,path,subset):
        merge_type = self.__get_subset_merge_type__(path)
        if merge_type == "pysam":
            try:
                file_object = pysam.AlignmentFile(path,"rb")
            except IOError:
                print "fail on path %s" % (path,)
                return
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
        for subset,file_obj in self._out_file_objects.items():
            file_obj.close()
        del self._out_file_objects

    def __handler_exit__(self,**kwargs):
        self.__close_all_out_files__()
        
#...happily ever after
