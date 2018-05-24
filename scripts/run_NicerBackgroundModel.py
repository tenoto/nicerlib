#!/usr/bin/env python

import sys 
import argparse
import subprocess

from nicerlib import nibgd 

if __name__ == '__main__':
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	help_message = """
Run the NICER Background Modeling 
(example) %s obsid_path outdir 
""" % (sys.argv[0])
	parser = argparse.ArgumentParser(
		description=help_message)
	parser.add_argument('obsid_path', metavar='obsid_path',type=str,	
		help='input obsid_path or @input_obsid_path.lst for multiple obsids.')
	parser.add_argument('outdir',metavar='outdir',type=str,
		help='output directory path.')
	parser.add_argument('--exposure_threshold',dest='exposure_threshold',action='store',
		help='GTI Exposure threshold (default 60 sec)', default=60.0)	
	parser.add_argument('--rmffile',dest='rmffile',action='store',
		help='rmffile', default="/Users/enoto/work/niresp/nicer_v1.02.rmf")		
	parser.add_argument('--arffile',dest='arffile',action='store',
		help='arffile', default="/Users/enoto/work/niresp/ni_xrcall_onaxis_v1.02.arf")			
#	parser.add_argument('--title',dest='title',action='store',
#		help='plot title', default=None)	
#	parser.add_argument('--outbase',dest='outbase',action='store',
#		help='outbase', default=None)	
	args = parser.parse_args()
	print(args)

	nibgd = nibgd.NicerBackgroundModel(args.obsid_path,args.outdir,
		args.exposure_threshold,args.rmffile,args.arffile)
	nibgd.run()

