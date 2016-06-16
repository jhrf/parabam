import pysam
import parabam
import time
import sys
import os
import gc
import shutil

from itertools import izip
from multiprocessing import Queue

class SubsetCore(object):

    def __init__(self,object constants,**kwargs):
        self.constants = constants

    def __pre_run_routine__(self,iterator,**kwargs):
        super(SubsetCore,self).__pre_run_routine__(iterator)
        self._subsets = self.constants.subsets
        self._counts = {}
        self._temp_paths = {}
        self._temp_objects = {}

        for subset in self._subsets:
            self._temp_paths[subset] = self.__get_temp_path__(subset)
            self._temp_objects[subset] = \
                    self.__get_temp_object__(self._temp_paths[subset])
            self._counts[subset] = 0

    def __get_results__(self):
        results = {}
        results["temp_paths"] = self._temp_paths
        results["counts"] = self._counts
        results["system"] = self._system

        return results

    def __write_to_subset_bam__(self,subset_type,read):
        self._counts[subset_type] += 1
        self._temp_objects[subset_type].write(read)

    def __post_run_routine__(self,**kwargs):
        super(SubsetCore,self).__post_run_routine__()
        for subset,fileobj in self._temp_objects.items():
            fileobj.close()


class Task(SubsetCore,parabam.command.Task):

    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants,**kwargs):
        
        parabam.command.Task.__init__(self,parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    statusqu=statusqu,
                                    task_size=task_size,
                                    constants=constants)
        SubsetCore.__init__(self,constants)

    def __handle_rule_output__(self,rule_output,read):             
        if type(rule_output) == bool:
            if rule_output:
                self.__write_to_subset_bam__(self._constants.subsets[0],read)          
        elif type(rule_output) == list:
            for subset,cur_read in rule_output:
                self.__write_to_subset_bam__(subset,cur_read)
        elif type(rule_output) == None:
            pass

class PairTask(SubsetCore,parabam.command.PairTask):
    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants):
        
        parabam.command.PairTask.__init__(self,parent_bam=parent_bam,
                                                inqu=inqu,
                                                outqu=outqu,
                                                statusqu=statusqu,
                                                task_size=task_size,
                                                constants=constants)
        SubsetCore.__init__(self,constants)

    def __handle_rule_output__(self,rule_output,read):
        for subset,cur_read in rule_output:
            self.__write_to_subset_bam__(subset,cur_read)

class ByCoordTask(SubsetCore,parabam.command.ByCoordTask):

    def __init__(self,parent_bam,inqu,outqu,statusqu,task_size,constants,**kwargs):
        
        parabam.command.ByCoordTask.__init__(self,parent_bam=parent_bam,
                                    inqu=inqu,
                                    outqu=outqu,
                                    statusqu=statusqu,
                                    task_size=task_size,
                                    constants=constants)
        SubsetCore.__init__(self,constants)

    def __handle_rule_output__(self,rule_output,read):     
        write_to_subset = self.__write_to_subset_bam__        
        for read in rule_output:
            write_to_subset(self._constants.subsets[0],read)

class Handler(parabam.command.Handler):

    def __init__(self,object parent_bam, object output_paths,object inqu,
                object constants,object pause_qus,dict out_qu_dict):
        
        super(Handler,self).__init__(parent_bam = parent_bam,
                                        output_paths = output_paths,
                                        inqu=inqu,
                                        constants=constants,
                                        pause_qus=pause_qus,
                                        out_qu_dict=out_qu_dict)

        self._subsets = constants.subsets
        for subset in self._subsets:
            self._stage_stores[subset] = []

    def __new_package_action__(self,new_package):
        super(Handler,self).__new_package_action__(new_package)
        results = new_package.results
        for subset in self._subsets:
            self._stage_stores[subset].append(\
                                (results["counts"][subset],
                                 results["temp_paths"][subset],
                                 new_package.sequence_id,))

    def __periodic_action__(self,iterations):
        super(Handler,self).__periodic_action__(iterations)

        for subset in self._subsets:
            if self.__test_stage_store__(subset):
                self.__add_merge_task__(results=self._stage_stores[subset],subset_type=subset)
                self.__update_stage_store__(subset)

    def __handler_exit__(self, **kwargs):
        super(Handler,self).__handler_exit__()
        #Kill the handlers handler

        for subset in self._subsets:
            #merge task 
            self.__add_merge_task__(self._stage_stores[subset],subset)
        self.__write_counts_csv__()
    
    def __add_merge_task__(self,results,subset_type):
        res = parabam.merger.MergePackage(results=results,
                                          subset_type=subset_type,)
        self._out_qu_dict["merge"].put(res)

    def __write_counts_csv__(self):
        if self._constants.output_counts:
            count_path = "./subset_counts_csv_%d.csv" % (time.time(),)
            self.__standard_output__(\
                    "- Outputting counts to following file: %s" % (count_path,))

            with open(count_path, "w") as count_file:
                key_order,header = self.__get_count_header__()
                count_file.write(header)
                for line in self.__generate_count_line__(key_order):
                    count_file.write(line)

    def __generate_count_line__(self,key_order):
        for source,stats in self._stats.items():
            line = ""
            line += "%s" % (source,)
            for key in key_order:
                line += ",%s" % (stats[key],)
            line += "\n"
            yield line

    def __get_count_header__(self):
        header = "Sample"
        keys = []
        for source, stats in self._stats.items():
            for key, val in stats.items():
                header += ","
                header += key
                keys.append(key)
            break
        header += "\n"
        return keys, header

class Subset(parabam.command.Interface):
    """The interface to parabam subset.
    Users will primarily make use of the ``run`` function."""

    def __init__(self, **kwargs):
        super(Subset,self).__init__(instance_name = "parabam subset", **kwargs)

    def run_cmd(self):
        module,rule,constants =\
                      self.__get_module_and_vitals__(self.cmd_args.rule)

        if hasattr(module,"get_subset_types"):
            subsets = module.get_subset_types()
        else:
            subsets = ["subset"]

        self.run(
            input_paths= self.cmd_args.input,
            constants = constants,
            rule = rule,
            subsets= subsets,
            fetch_region = self.cmd_args.region,
            output_counts= self.cmd_args.counts)
    
    def run(self,input_paths,
            constants,
            rule,
            subsets,
            fetch_region=None,
            output_counts=False,
            **kwargs):

        ''' Docstring! '''
        args = dict(locals())
        del args["self"]
        results = super(Subset,self).run(**args)
        return results

    def __setup_cmd_line_run__(self):

        super(Subset,self).__setup_cmd_line_run__()
        self.ensure_unique_output = self.cmd_args.u

    def __get_queue_names__(self,**kwargs):
        queues = ["merge","main"]
        return queues

    def __get_destroy_handler_order__(self):
        return [Handler,parabam.merger.Handler]

    def __get_handler_bundle__(self,**kwargs):
        handler_bundle = {}

        #queues transformed by leviathon
        handler_bundle[Handler] = {"inqu":"main", 
                                   "out_qu_dict":["merge"]}

        handler_bundle[parabam.merger.Handler] = {"inqu":"merge",
                                                  "out_qu_dict":[]}
        return handler_bundle

    def __instalise_final_output__(self,**kwargs):
        return {}

    def __get_output_paths__(self,
                             input_path,
                             final_output_paths,
                             subsets,
                             **kwargs):
        
        output_paths = {input_path:{}}
        for salt,subset in enumerate(subsets):
            output_paths[input_path][subset] =\
                                     self.__get_path__(input_path,
                                                         subset,
                                                         salt)
        return output_paths

    def __get_path__(self,input_path,subset,salt):
        head,ext = os.path.splitext(input_path)
        head,tail = os.path.split(head)
        unique = ""

        if self.ensure_unique_output:
            unique = "_%d%d" % (time.time(),salt)

        return os.path.join(".",self.temp_dir, "%s_%s%s%s" \
                                    % (tail,subset,unique,ext))
            
    def __get_task_class__(self,**kwargs):
        if self.pair_process:
            return PairTask
        elif self.coord_process:
            return ByCoordTask
        else:
            return Task

    def get_parser(self):
        parser = self.default_parser()

        parser.add_argument('-u',action="store_true",default=False,
            help="The files will be appended with a unique identifier."\
                  " In the format <input_name>_<unique_id>_<subset_type>.bam")
        parser.add_argument('--counts',action="store_true",default=False,
            help="The amount of reads in each subset will be output as"\
                    " a .CSV alongside the .BAM file.")

        return parser 

#...happily ever after
