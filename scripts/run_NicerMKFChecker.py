#!/usr/bin/env python

import sys 
import argparse
import subprocess

from nicerlib import nimkf 


if __name__ == '__main__':
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	help_message = """
Run the basic NICER mkf checker.
(example) %s mkf file
""" % (sys.argv[0])
	parser = argparse.ArgumentParser(
		description=help_message)
	parser.add_argument('mkffile', metavar='mkffile',type=str,	
		help='input mkffile')
	parser.add_argument('outdir', metavar='outdir',type=str,	
		help='output directory')	
	args = parser.parse_args()
	print(args)

	nimkf = nimkf.NicerMKF(args.mkffile, args.outdir)
	nimkf.run()

