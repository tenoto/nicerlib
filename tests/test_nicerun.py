#!/usr/bin/env python

import argparse

import sys 
sys.path.append('../nicerlib/')
import nicerun

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run the basic NICER analysis process.')
	parser.add_argument('obsid_path', metavar='obsid_path',type=str,	
		help='input obsid path')
#	parser.add_argument('--outdir', dest='outdir_main', 
#		type=str, default=__outdir__,
#		help='output main directory name.')		
	args = parser.parse_args()

	niobs = nicerun.NicerObservation(args.obsid_path)
	niobs.show_parameter()