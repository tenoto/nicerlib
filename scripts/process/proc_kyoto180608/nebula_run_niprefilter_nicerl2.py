#!/usr/bin/env python

import os 
import sys 
import yaml 
import argparse

fparam_default = '%s/scripts/process/proc_kyoto180608/proc_param.yaml' % os.getenv('NICER_SOFT_PATH')

help_message = """
(example) %s --indir indir --fparam fparam
""" % sys.argv[0]
parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('obsid', metavar='obsid',type=str,help='ObsID') 
parser.add_argument('outdir', metavar='outdir',type=str,
	help='output directory.') 
parser.add_argument('--fparam', metavar='fparam',type=str,        
	help='yaml file for input parameters.',
	default=fparam_default)
parser.add_argument("-r", "--recreate",action="store_true",
	dest="flag_recreate", default=False,
	help='recreate flag.')
args = parser.parse_args()
print(args)

obsid_path = '/Users/enoto/work/drbv1/reporitory/heasarc/data/nicer/data/obs/*/%s' % args.obsid

if not os.path.exists(args.outdir):
	cmd = 'mkdir -p %s' % args.outdir
	print(cmd);os.system(cmd)

proc_outdir = '%s/%s' % (args.outdir,args.obsid)
print(proc_outdir)

if os.path.exists(proc_outdir) and not args.flag_recreate:
	print("file %s has already existed." % proc_outdir)
	print("............skipped.")
	exit()

if args.flag_recreate:
	cmd = 'rm -rf %s' % proc_outdir
	print(cmd);os.system(cmd)

cmd  = 'cp -r %s %s' % (obsid_path,args.outdir)
print(cmd);os.system(cmd)

cmd  = 'run_niprefilter_nicerl2.py '
cmd += '--indir %s ' % proc_outdir
cmd += '--fparam %s ' % args.fparam
print(cmd);os.system(cmd)