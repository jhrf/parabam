#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os,argparse
import pysam
import time
import parabam
import gzip

import re

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
                      object constants,object pause_qus,dict out_qu_dict):

        super(Handler,self).__init__(parent_bam,output_paths,inqu,
                                        constants,pause_qus,
                                        out_qu_dict,report=False)
        
        self._user_subsets = list(constants.user_subsets)

        self._eof_signature = '\x1F\x8B\x08\x04\x00\x00\x00\x00\x00\xFF'+\
                              '\x06\x00\x42\x43\x02\x00\x1B\x00\x03\x00'+\
                              '\x00\x00\x00\x00\x00\x00\x00\x00'
        self._header_signature = '\x1F\x8B\x08\x04\x00\x00\x00\x00\x00\xFF'

        self._out_file_objects = self.__get_out_file_objects__()
        self._merged = 0

    def __get_out_file_objects__(self):
        file_objects = {}
        output_paths = self._output_paths

        file_objects = {}
        for subset in self._user_subsets:
            output_path = output_paths[self._parent_bam.filename][subset]
            cur_file_obj = open(output_path,"wb")
            self.__write_header_to_file__(cur_file_obj)
            file_objects[subset] = cur_file_obj

        return file_objects

    def __write_header_to_file__(self,cur_file_obj):
        header_location = self.__get_header_location__\
                                    (self._parent_bam.filename)

        with open(self._parent_bam.filename,"rb") as header_file:
            cur_file_obj.write(header_file.read(header_location))
            cur_file_obj.flush()

    def __get_header_location__(self,header_path):
        header_signature = self._header_signature
        i = 0

        with open(header_path,"rb") as header_file:
            read_bytes = ""

            while True:
                read_bytes += header_file.read(2048)
                header_count = read_bytes.count(header_signature)

                if header_count > 1:
                    second_header = [i.start() for i in \
                                re.finditer(header_signature,read_bytes)][1]
                    return second_header

                i +=1
                if i > 50:
                     # TODO: Handle this error better
                    print "Header not found."
                    print "Probably not a BAM file"
                    sys.exit(0)

    def __periodic_action__(self,iterations):
        if self._destroy:
            try:
                pack = self._inqu.get(True,10)
                self._inqu.put(pack)
            except Queue2.Empty:
                self._finished = True

    def __new_package_action__(self,new_package,**kwargs):
        #NewPackge is of type parabam.merger.MergePackage
        subset_type = new_package.subset_type
        
        for merge_count,merge_path in new_package.results:
            if merge_count > 0:
                with open(merge_path,"rb") as merge_file:
                    merge_file.seek(self.__get_header_location__(merge_path))

                    for binary_data in \
                            self.__get_binary_from_file__(merge_file):

                        self._out_file_objects[subset_type].write(binary_data)
                        self._out_file_objects[subset_type].flush()

            #time.sleep(20)
            os.remove(merge_path)
            self._merged += 1
        time.sleep(.5)

    def __get_binary_from_file__(self,open_file):
        prev_data = ""
        i = 0

        while True:
            binary_data = open_file.read(2**19)
            if not i == 0:
                combined = prev_data+binary_data
                if combined[len(combined)-28:] == self._eof_signature:
                    # print "Last"
                    # print combined[len(combined)-1024:].encode("hex")
                    yield combined[:len(combined)-28]
                    return
                else:
                    #if i == 1:
                        # print "First"
                        # print prev_data[:1024].encode("hex")
                    yield prev_data

            prev_data = binary_data
            i += 1

    def __close_all_out_files__(self):
        for subset,file_obj in self._out_file_objects.items():
            file_obj.write(self._eof_signature)
            file_obj.flush()
            file_obj.close()
        del self._out_file_objects

    def __handler_exit__(self,**kwargs):
        self.__close_all_out_files__()
        
#...happily ever after
