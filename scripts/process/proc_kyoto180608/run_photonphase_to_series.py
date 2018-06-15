#!/usr/bin/env python

import os 
import sys 
import yaml 
import glob
import argparse

fparam_default = '%s/scripts/process/proc_kyoto180608/proc_param.yaml' % os.getenv('NICER_SOFT_PATH')

help_message = """
(example) %s --indir indir --fparam fparam
""" % sys.argv[0]
parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('obsid_base', metavar='obsid',type=str,help='ObsID') 
parser.add_argument('outdir', metavar='outdir',type=str,
	help='output directory.') 
parser.add_argument('parfile', metavar='parfile',type=str,        
	help='parfile for photonphase.')
parser.add_argument('--fparam', metavar='fparam',type=str,        
	help='yaml file for input parameters.',
	default=fparam_default)
parser.add_argument('--ephem', metavar='ephem',type=str,        
	help='ephemeris (e.g., DE421) for photonphase.',
	default="")

args = parser.parse_args()
print(args)

param = yaml.load(open(args.fparam))
print(param['gtiexpr_basestr'])


def add_photonphase(obsid_path):
	obsid = os.path.basename(obsid_path).strip()
	print("---------------------------------------------")
	print("obsid:%s" % obsid)
	orbfile  = '%s/auxil/ni%s.orb.gz' % (obsid_path,obsid)
	clsrcevt = '%s/proc/speclc/prod/ni%s_%s_clscr_gtisel.evt' % (obsid_path,obsid,param['gtiexpr_basestr'])

	if not os.path.exists(clsrcevt):
		print("file %s does not exist. skipped." % clsrcevt)
		return -1

	photonphase_evt = clsrcevt.replace('.evt','_brph.evt')	
	photonphase_log = clsrcevt.replace('.evt','_brph.log')		
	photonphase_png = clsrcevt.replace('.evt','_brph.png')			

	cmd = 'cp %s %s' % (clsrcevt,photonphase_evt)
	print(cmd);os.system(cmd)	

	cmd  = 'photonphase --addphase '
	cmd += '--ephem %s ' % args.ephem
	cmd += '--orbfile %s ' % orbfile 
	cmd += '--plotfile %s ' % photonphase_png
	cmd += '%s %s >& ' % (photonphase_evt,args.parfile)
	cmd += '%s ' % photonphase_log
	print(cmd);os.system(cmd)	

	return 0 

for obsid_path in glob.glob('%s/%s*' % (args.outdir,args.obsid_base)):
	add_photonphase(obsid_path)
	
