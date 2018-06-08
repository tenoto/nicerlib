#!/usr/bin/env python

__name__    = 'run_niprefilter_nicerl2'
__author__  = 'Teru Enoto'
__version__ = 'proc_kyoto180608'
__date__    = '2018 June 8'

"""
run following commands (similar to nicerl2)
1. niprefilter2 - create augmented NICER-specific filter file
2. nicercal - apply standard NICER calibration
3. nimaketime - create standard screening good time intervals
4. nicermergeclean - combine per-MPU data and filter/screen
"""

import os 
import sys 
import yaml 
import argparse

# ==============================
# Get input parameters 
# ==============================
help_message = """
(example) %s obsid_path.lst 
""" % sys.argv[0]

parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('--indir', metavar='indir',type=str,        
	help='Input directory name. The directory should be a single NICER observation directory, which in turn contains xti/{events_uf,events_cl,hk,auxil} subdirectories.')
parser.add_argument('--fparam', metavar='fparam',type=str,        
	help='yaml file for input parameters.')
parser.add_argument("-r", "--recreate",action="store_true",
	dest="flag_recreate", default=False,
	help='recreate flag.')
args = parser.parse_args()
print(args)

if not os.path.exists(args.fparam):
	sys.stderr.write('error: file %s does not exist.' % args.fparam)
	quit()

param = yaml.load(open(args.fparam))
print(param)

# ==============================
# Prepare
# ==============================
obsid = os.path.basename(args.indir)
fname_mkffile = '%s/auxil/ni%s.mkf.gz' % (args.indir,obsid)

dir_proc = '%s/proc' % args.indir
dir_log  = '%s/log' % dir_proc
cmd = 'mkdir -p %s' % dir_log
print(cmd);os.system(cmd)

fname_gtifile = '%s/ni%s_%s.gti' % (dir_proc,obsid,param['gtiexpr_basestr'])

# ==============================
# 1. niprefilter2 - create augmented NICER-specific filter file
# ==============================
cmd = 'gunzip %s' % fname_mkffile
print(cmd);os.system(cmd)

fname_log_niprefilter2 = '%s/1_niprefilter2.log' % dir_log
cmd = 'niprefilter2 clobber=yes indir=%s infile=%s outfile=INFILE >& %s' % (args.indir,fname_mkffile.replace('.gz',''),fname_log_niprefilter2)
print(cmd);os.system(cmd)

cmd = 'gzip %s' % fname_mkffile.replace('.gz','')
print(cmd);os.system(cmd)

# ==============================
# 2. nicercal - apply standard NICER calibration
# ==============================
fname_log_nicercal = '%s/2_nicercal.log' % dir_log
fname_outfilefile = '%s/ufalist.lis' % dir_proc
cmd = 'nicercal indir=%s outfilefile=%s >& %s' % (args.indir,fname_outfilefile,fname_log_nicercal)
print(cmd);os.system(cmd)

cmd = 'gzip %s' % fname_mkffile.replace('.gz','')
print(cmd);os.system(cmd)

# ==============================
# 3. nimaketime - create standard screening good time intervals
# ==============================
fname_log_nimaketime = '%s/3_nimaketime.log' % dir_log
cmd = 'nimaketime infile=%s outfile=%s expr="%s" outexprfile="NONE" >& %s' % (fname_mkffile, fname_gtifile, param['nimaketime_expr'], fname_log_nimaketime)
print(cmd);os.system(cmd)
exit()

# ==============================
# 4. nicermergeclean - combine per-MPU data and filter/screen
# ==============================
fname_ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (args.indir,obsid)
fname_clevt  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (args.indir,obsid)
cmd  = 'nicermergeclean infiles=@%s ' % fname_outfilefile
cmd += 'ufafile=%s ' % fname_ufaevt
cmd += 'clfile=%s ' % fname_clevt
cmd += 'gtifile=%s ' % fname_gtifile
print(cmd);os.system(cmd)











