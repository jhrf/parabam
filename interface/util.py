import argparse
import os
import shutil
import pysam
import time
def default_parser():
	parser = argparse.ArgumentParser(conflict_handler='resolve')

	parser.add_argument('-p',type=int,nargs='?',default=8
		,help="The amount of processors you wish the program to use."\
				" Ignored if argument 's' is provided")
	parser.add_argument('-c',type=int,nargs='?',default=25000
		,help="How many arguments the preprocesor should store before"\
		"sending to the analysis module")
	parser.add_argument('-v',action="store_true",default=False
		,help='If set the program will output more information to terminal')
	parser.add_argument('--instruc','-i',metavar='INSTRUCTION',required=True
		,help='The instruction file, written in python, that we wish'\
		'to carry out on the input BAM.')
	parser.add_argument('--input','-b',metavar='INPUT', nargs='+',required=True
		,help='The file(s) we wish to operate on. Multipe entries should be separated by a single space')
	parser.add_argument('--output','-o',metavar='OUTPUT', nargs='+',required=False
		,help='The name of the output that we wish to create. Must be same amount of space \
		separated entries as INPUT.')

	return parser

def get_bam_basename(path):
	base = os.path.basename(path)
	if "." in base:
		base = base.rpartition(".")[0]
	return base

def create_dirs(dirs):
	for i in dirs:
		if i == "": continue
		os.makedirs(i)

def die_gracefully(tempDir):
	if not tempDir == "." or not tempDir == "": 
		shutil.rmtree(tempDir)

def get_unique_temp(filID,tempDir="."):
	return "%s/%sTEMP%d.bam" % (tempDir,filID,int(time.time()),)

def create_subset_bam(master_file_path,subset_file_path):
	master = pysam.Samfile(master_file_path,"rb")
	subset = pysam.Samfile(subset_file_path,"wb",template=master)

	master.close()
	subset.close()