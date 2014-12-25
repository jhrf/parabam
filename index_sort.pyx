import pysam
import parabam
import pdb
import time
import sys
import os
import gc
import shutil
import gzip

from multiprocessing import Queue
from parabam.subset import PairTaskSubset

class TaskIndex(PairTaskSubset):

    def __init__(self,object task_set,object outqu,object curproc,object destroy,object const,str source):
        super(TaskIndex, self).__init__(task_set=task_set,
                                        outqu=outqu,
                                        curproc=curproc*len(const.subset_types),
                                        destroy=destroy,
                                        const=const,
                                        source=source)
        
    def __handle_task_set__(self,task_set,master):
        engine = self.__engine__
        user_constants = self.const.user_constants
        subset_types = self.const.subset_types

        subset_write = self.__write_to_subset_bam__

        loners = {}
        pairs = {}

        #Find pairs in index chunk
        for index in task_set:
            master.seek(int(index))
            read = master.next()
            read1,read2 = self.__query_loners__(read,index,loners)
            if read1:
                pairs = {read1.qname:(read1,read2)}
        
        #Assign reads that are still loners to loners file
        for read,index in loners.items():
            subset_write("index",read)
        del loners

        #Send pairs onto the subset decision
        task_set

class ProcessorIndex(parabam.core.Processor):

    def __init__(self,object outqu,object const,object TaskClass,object debug=False):
        # if const.fetch_region:
        #   debug = False 
        super(ProcessorIndex,self).__init__(outqu,const,TaskClass,task_args,debug)

    def __get_master_bam__(self,master_file_path):
        return gzip.open(,"rb")

    def __add_to_collection__(self,master,alig,collection):
        collection.append(alig)

    def __get_next_alig__(self,master_bam):
        for alig in master_bam:
            yield alig

    def __pre_processor__(self,master_file_path):
        pass

class IndexInterface(parabam.core.Interface):
    """The interface to parabam subset.
    Users will primarily make use of the ``run`` function."""

    def __init__(self,temp_dir,exe_dir):
        super(Interface,self).__init__(temp_dir,exe_dir)

    def run_cmd(self,parser):
        print "[Error] This module cannot be run from the command line."
        pass
    
    def run(self,master_file_paths,index_path,proc,
            chunk,const,user_constants,user_engine,
            subset_types,task_class,
            source,verbose=False,debug=False,):

        run = 0
        while True:
            basename,ext = os.path.splitext(os.path.basename(master_file_path[source]))
            output_paths = self.__get_output_paths__(source,subset_types,basename,run)
            subset_types.append("index")
                    
            const.add("task_class",task_class)
            const.add("index_path",index_path)
            const.add("source",source)

            task_qu = Queue()

            # outqu,const,TaskClass,debug=False
            processors = [ProcessorIndex(Queue,const,)]
            handlers = self.__create_handlers__(task_qu,const)

            if verbose == 1: 
                update_interval = 199
            else:
                update_interval = 1

            lev = parabam.core.Leviathon(processors,handlers,update_interval)
            lev.run()
            del lev
            
            gc.collect()
            run += 1

            pdb.set_trace()

        return final_files

    def __get_output_paths__(self,source,subset_types,basename,run):
        paths = {source:{}}
        for subset in subset_types:
            paths[source][subset] = "%s/index_%s_%s_%d.bam" % (self._temp_dir,subset,basename,run)
        paths[source]["index"] = "%s/leftovers_%s_%s_%d.gz" % (self._temp_dir,subset,basename,run)

        return paths

#...happily ever after
