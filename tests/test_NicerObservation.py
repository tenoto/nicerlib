#!/usr/bin/env python

import argparse
import subprocess

from nicerlib import nicerun 

if __name__ == '__main__':
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	parser = argparse.ArgumentParser(description='Run the basic NICER analysis process (NicerObservation).')
	parser.add_argument('obsid_path', metavar='obsid_path',type=str,	
		help='input obsid path')
	args = parser.parse_args()

	niobs = nicerun.NicerObservation(args.obsid_path)
	niobs.run_nicerl2()
	niobs.extract_overonly_event()
