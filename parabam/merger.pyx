#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os
import time
import re
import parabam
import pysam

import Queue as Queue2

from collections import namedtuple
from itertools import izip

from pprint import pprint as ppr


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
        
        self._subsets = list(constants.subsets)

        #CONSTANTS
        self._CLUMP_THRESH = 1000
        self._EOF_SIGNATURE = '\x1F\x8B\x08\x04\x00\x00\x00\x00\x00\xFF'+\
                              '\x06\x00\x42\x43\x02\x00\x1B\x00\x03\x00'+\
                              '\x00\x00\x00\x00\x00\x00\x00\x00'
        self._Clump = namedtuple("Clump","sequence_ids merge_tuples")                              

        #VARIABLES

        self._out_file_objects = self.__get_out_file_objects__()
        self._subset_has_header = \
                self.__get_subset_dict__(payload=(lambda:False))
        self._merged = 0

        self._preserve_sequence = not self._constants.pair_process
        self._waiting_files = self.__get_subset_dict__(payload=dict)
        self._merge_stage = self.__get_subset_dict__()
        self._previous_merge = -1
        self._out_of_sequence_package = -2

        self._total_clumps = 0

        self._header_path = self.__get_header_sam_path__()

    def __get_header_sam_path__(self):
        header = pysam.view(*["-H",self._parent_bam.filename])
        time_of_creation = int(time.time())
        header_path = os.path.join(self._constants.temp_dir,
                                   "header_%d.bam" % (time_of_creation,))
        with open(header_path,"w") as header_sam_file:
            header_sam_file.write("".join(header))
            header_sam_file.flush()
        return header_path

    def __get_subset_dict__(self,payload=list):
        return dict( ( (subset,payload(),) for subset in self._subsets) )

    def __get_out_file_objects__(self):
        file_objects = {}
        output_paths = self._output_paths

        file_objects = {}
        for subset in self._subsets:
            output_path = output_paths[self._parent_bam.filename][subset]
            file_objects[subset] = self.__get_bam_file_obj__(output_path)

        return file_objects

    def __get_bam_file_obj__(self,output_path):
        return open(output_path,"wb")

    def __get_header_location__(self,header_path):
        header_bytes = pysam.view(*["-Hb",header_path])
        header_str = "".join(header_bytes).encode("hex")
        return (len(header_str) / 2) - len(self._EOF_SIGNATURE)

    def __periodic_action__(self,iterations):
        if self._destroy:
            try:
                pack = self._inqu.get(True,10)
                self._inqu.put(pack)
            except Queue2.Empty:
                self._finished = True

    def __new_package_action__(self,new_package,**kwargs):
        #new_package is of type parabam.merger.MergePackage
        subset = new_package.subset_type
        self.__add_to_waiting_files__(subset,new_package.results)
        self.__process_waiting_files__(subset)

    def __add_to_waiting_files__(self,subset,results):
        for merge_count,merge_path,sequence_id in results:
            merge_id = sequence_id
            if merge_id == -1:
                merge_id = self._out_of_sequence_package
                self._out_of_sequence_package += (-1)

            self._merge_stage[subset].append(merge_id)
            self._waiting_files[subset][merge_id] = (merge_count,merge_path)
        self._merge_stage[subset].sort()

    def __process_waiting_files__(self,subset,force=False):
        #TODO: Better name for valid_sequence_ids
        valid_sequence_ids = self.__get_valid_sequence_ids__(subset) 
        if len(valid_sequence_ids) > 0:

            merge_tuples = [ self._waiting_files[subset][sequence_id]\
                             for sequence_id in valid_sequence_ids ]

            clumps,clumped_sequence_ids = \
                    self.__clump_merge_files__(merge_tuples,
                                               valid_sequence_ids,
                                               force)

            if len(clumps) > 0:
                for merge_count,merge_path in clumps:
                    self.__dump_to_BAM_file__(merge_path,subset)
            
                for sequence_id in clumped_sequence_ids:
                    del self._waiting_files[subset][sequence_id]
                    self._merge_stage[subset].remove(sequence_id)

            del clumps
            del valid_sequence_ids
            del clumped_sequence_ids

    def __clump_merge_files__(self,merge_tuples,sequence_ids,force=False):
        all_above_thresh = True
        any_above_thresh = False
        total_reads = 0

        for merge_count,merge_path in merge_tuples:
            total_reads += merge_count

            if merge_count >= self._CLUMP_THRESH:
                any_above_thresh = True
            else:
                all_above_thresh = False

        if all_above_thresh:
            #no clump required
            return merge_tuples,list(sequence_ids)
        elif (total_reads < self._CLUMP_THRESH and not force) \
            or (force and total_reads == 0):
            #clump not possible
            return [],[]
        elif not any_above_thresh or force:
            #one large clump
            return [ self.__get_clump__(merge_tuples) ],list(sequence_ids)
        else:
            #several clumps
            return self.__create_clumps__(merge_tuples,list(sequence_ids))
    
    def __create_clumps__(self,merge_tuples,sequence_ids):
        clumped_sequence_ids = []
        new_merge_tuple = []

        current_clump = self.__new_clump__()
        count = 0

        for (merge_count,merge_path),sequence_id in \
                                            izip(merge_tuples,sequence_ids):
            if count > self._CLUMP_THRESH:
                new_merge_tuple.append(self.__get_clump__(\
                                            current_clump.merge_tuples))
                clumped_sequence_ids.extend(current_clump.sequence_ids)
                current_clump = self.__new_clump__()
                count = 0

            current_clump.merge_tuples.append((merge_count,merge_path))
            current_clump.sequence_ids.append(sequence_id)
            count += merge_count

        return new_merge_tuple,clumped_sequence_ids

    def __new_clump__(self):
        return self._Clump(sequence_ids=[],merge_tuples=[])

    def __get_clump__(self,merge_tuples):

        # This function may seem overly complicated but it is 
        # neccessary to handle input to pysam.cat in this way.
        #
        # The function protects pysam.cat from catting too many
        # files at once and from catting only a single file
        if len(merge_tuples) == 1:
            return merge_tuples[0]
        else:

            def get_tuple_generator(merge_tuples):
                for count,path in merge_tuples:
                    yield count,path
                return

            new_tuples = []
            total_clumped = 0
            merge_tuples_generator = get_tuple_generator(merge_tuples)

            while total_clumped < len(merge_tuples):

                clump_path = self.__get_clump_path__()
                clump_count = 0
                
                remove_paths = []
                cat_paths = []

                for count,path in merge_tuples_generator:
                    total_clumped += 1
                    remove_paths.append(path)
                    if count > 0:
                        cat_paths.append(path)
                        clump_count += count

                    if len(cat_paths) >= 250:
                        break
                

                if len(cat_paths) > 1:
                    cat_paths.insert(0,"-o%s" % (clump_path,))
                    cat_paths.insert(0,"-h%s" % (self._header_path,))
                    pysam.cat(*cat_paths)
                    new_tuples.append((clump_count,clump_path,))

                elif len(cat_paths) == 1:
                    # Rare case where a sub-clump was comprised of empty BAM 
                    # files apart from one BAM file with reads

                    new_tuples.append((clump_count,cat_paths[0]),)
                    remove_paths.remove(cat_paths[0])
                

                for path in remove_paths:
                    os.remove(path)

            return self.__get_clump__(new_tuples)

    def __get_clump_path__(self):
        clump_path = os.path.join(self._constants.temp_dir,
                              "merge_temp_%d-%d.bam" % (self._total_clumps,
                                                         time.time(),))
        self._total_clumps += 1
        return clump_path

    def __get_valid_sequence_ids__(self,subset):
        if not self._preserve_sequence:
            return self._merge_stage[subset]
        else:
            ready_to_merge = []
            prev = self._previous_merge

            self._merge_stage[subset].sort()
            for sequence_id in self._merge_stage[subset]:
                
                if sequence_id-prev > 1:
                    break
                else:
                    ready_to_merge.append(sequence_id)
                    prev = sequence_id

            if len(ready_to_merge) > 0:
                self._previous_merge = ready_to_merge[-1]

            return ready_to_merge

    def __dump_to_BAM_file__(self,merge_path,subset):
        with open(merge_path,"rb") as merge_file:
            if self._subset_has_header[subset]:
                merge_file.seek(self.__get_header_location__(merge_path))
            else:
                self._subset_has_header[subset] = True

            for binary_data in \
                self.__get_binary_from_file__(merge_file):

                self._out_file_objects[subset].write(binary_data)
                self._out_file_objects[subset].flush()
        os.remove(merge_path)

    def __get_binary_from_file__(self,open_file):
        prev_data = ""
        i = 0

        while True:
            binary_data = open_file.read(2**19)
            if not i == 0:
                combined = prev_data+binary_data
                if combined[len(combined)-28:] == self._EOF_SIGNATURE:
                    yield combined[:len(combined)-28]
                    return
                else:
                    yield prev_data

            prev_data = binary_data
            i += 1

    def __close_all_out_files__(self):
        for subset,file_obj in self._out_file_objects.items():
            file_obj.write(self._EOF_SIGNATURE)
            file_obj.flush()
            file_obj.close()
        del self._out_file_objects

    def __handler_exit__(self,**kwargs):
        for subset in self._subsets:
            self.__process_waiting_files__(subset,force=True)
        self.__close_all_out_files__()
        
#...happily ever after
