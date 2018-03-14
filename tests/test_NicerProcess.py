#!/usr/bin/env python

import sys 
import argparse
import subprocess

from nicerlib import nicerun 


if __name__ == '__main__':
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	help_message = """
Run the basic NICER analysis process (NicerRun).
(example1) %s $NICER_DATA_PATH/obs2validate/1/1012040107 test_NicerProcess_fparam.yaml outdir 
(example2) %s @test_NicerProcess_obsid_path.txt test_NicerProcess_fparam.yaml outdir 
""" % (sys.argv[0],sys.argv[0])
	parser = argparse.ArgumentParser(
		description=help_message)
	parser.add_argument('obsid_path', metavar='obsid_path',type=str,	
		help='input obsid_path or @input_obsid_path.lst for multiple obsids.')
	parser.add_argument('fparam',metavar='fparam',type=str,
		help='input parameter yaml file for the pipeline.')
	parser.add_argument('outdir',metavar='outdir',type=str,
		help='output directory path.')
	parser.add_argument('--title',dest='title',action='store',
		help='plot title', default=None)	
	parser.add_argument('--outbase',dest='outbase',action='store',
		help='outbase', default=None)	
	args = parser.parse_args()
	print(args)

	niproc = nicerun.NicerProcess(args)
	niproc.run()

