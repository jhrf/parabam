#!/usr/bin/env python
#Once upon a time...

import pdb,sys,os
import time
import re
import parabam
import pysam

import Queue as Queue2


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
        self._subset_has_header = \
                self.__get_subset_dict__(payload=(lambda:False))
        self._merged = 0

        self._preserve_order = not self._constants.pair_process
        self._waiting_files = self.__get_subset_dict__(payload=dict)
        self._merge_stage = self.__get_subset_dict__()
        self._previous_merge = -1

        self._total_clumps = 0

    def __get_subset_dict__(self,payload=list):
        return dict( ( (subset,payload(),) for subset in self._user_subsets) )

    def __get_out_file_objects__(self):
        file_objects = {}
        output_paths = self._output_paths

        file_objects = {}
        for subset in self._user_subsets:
            output_path = output_paths[self._parent_bam.filename][subset]
            file_objects[subset] = self.__get_bam_file_obj__(output_path)

        return file_objects

    def __get_bam_file_obj__(self,output_path):
        return open(output_path,"wb")

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
        #new_package is of type parabam.merger.MergePackage
        subset = new_package.subset_type
        self.__add_to_waiting_files__(subset,new_package.results)
        self.__process_waiting_files__(subset)

    def __add_to_waiting_files__(self,subset,results):
        for merge_count,merge_path,order in results:
            self._merge_stage[subset].append(order)
            self._waiting_files[subset][order] = (merge_count,merge_path)
        self._merge_stage[subset].sort()

    def __process_waiting_files__(self,subset,force=False):
        order_check = self.__order_check__(subset) #TODO: Better name for order_check
        if len(order_check) > 0:
            merge_tuples = [ self._waiting_files[subset][order] for order in order_check ]

            solve_clumps = self.__clump_merge_files__(merge_tuples,force)
            
            if len(solve_clumps) > 0:
                for merge_count,merge_path in solve_clumps:
                    self.__dump_to_BAM_file__(merge_path,subset)
               
                for order in order_check:
                    del self._waiting_files[subset][order]
                    self._merge_stage[subset].remove(order)

            del solve_clumps
            del order_check

    def __clump_merge_files__(self,merge_tuples,force=False):
        all_above_thresh = True
        any_above_thresh = False
        total_reads = 0

        merge_paths = []

        for merge_count,merge_path in merge_tuples:
            total_reads += merge_count

            if merge_count >= 512:
                any_above_thresh = True
            else:
                all_above_thresh = False

        if all_above_thresh:
            #no clump required
            return merge_tuples
        elif total_reads < 2000 and not force:
            #clump not possible
            return []
        elif not any_above_thresh:
            #one large clump
            return [ (total_reads,self.__get_clump__(merge_tuples),) ]
        else:
            #several clumps
            return self.__clump_merge_files__(self.__pair_small_files__(merge_tuples),force)
    
    def __pair_small_files__(self,merge_tuples):
        new_tuples = []
        for t1 in xrange(0,len(merge_tuples),2):
            if t1+1 < len(merge_tuples):
                pair1 = merge_tuples[t1]
                pair2 = merge_tuples[t1+1]

                if pair1[0] < 512 or pair2[0] < 512:
                    new_tuples.append( self.__get_clump__([pair1,pair2]) )
                else:
                    new_tuples.extend((pair1,pair2,))
            else:
                new_tuples.append(merge_tuples[t1])
        return new_tuples

    def __get_clump__(self,merge_tuples):
        clump_path = os.path.join(self._constants.temp_dir,"merge_temp_%d-%d.bam" \
                                                    % (self._total_clumps,time.time(),))
        clump_file = pysam.AlignmentFile(clump_path,"wb",header=self._parent_bam.header)
        clump_count = 0
                
        for count,path in merge_tuples:
            clump_count += count

            if count > 0:
                small_bam = pysam.AlignmentFile(path,"rb")
                for read in small_bam.fetch(until_eof=True):
                    clump_file.write(read)
                small_bam.close()
            os.remove(path)

        clump_file.close()
        self._total_clumps += 1
        return clump_count,clump_path

    def __order_check__(self,subset):
        if not self._preserve_order:
            return self._merge_stage[subset]
        else:
            ready_to_merge = []
            prev = self._previous_merge

            self._merge_stage[subset].sort()
            for order in self._merge_stage[subset]:
                
                if order-prev > 1:
                    break
                else:
                    ready_to_merge.append(order)
                    prev = order

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
                if combined[len(combined)-28:] == self._eof_signature:
                    yield combined[:len(combined)-28]
                    return
                else:
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
        for subset in self._user_subsets:
            self.__process_waiting_files__(subset,force=True)
        self.__close_all_out_files__()
        
#...happily ever after
