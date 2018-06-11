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
import astropy.io.fits as pyfits 

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
#parser.add_argument("-r", "--recreate",action="store_true",
#	dest="flag_recreate", default=False,
#	help='recreate flag.')
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
fname_outfilefile = '%s/2_nicercal_ufalist.lis' % dir_log
cmd = 'nicercal indir=%s outfilefile=%s >& %s' % (args.indir,fname_outfilefile,fname_log_nicercal)
print(cmd);os.system(cmd)

cmd = 'gzip %s' % fname_mkffile.replace('.gz','')
print(cmd);os.system(cmd)

# ==============================
# 3. nimaketime - create standard screening good time intervals
# ==============================
fname_log_nimaketime = '%s/3_nimaketime.log' % dir_log
fname_expr_nimaketime = '%s/3_nimaketime.expr' % dir_log
fname_gtifile = '%s/ni%s_%s.gti' % (dir_proc,obsid,param['gtiexpr_basestr'])
cmd = 'nimaketime infile=%s outfile=%s expr="%s" outexprfile="%s" >& %s' % (fname_mkffile, fname_gtifile, param['nimaketime_expr'], fname_expr_nimaketime, fname_log_nimaketime)
print(cmd);os.system(cmd)

# ==============================
# 4. nicermergeclean - combine per-MPU data and filter/screen
# ==============================
fname_log_nicermergeclean = '%s/4_nicermergeclean.log' % dir_log
fname_ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa_%s.evt' % (args.indir,obsid,param['gtiexpr_basestr'])
fname_clevt  = '%s/xti/event_cl/ni%s_0mpu7_cl_%s.evt' % (args.indir,obsid,param['gtiexpr_basestr'])
cmd  = 'nicermergeclean infiles=@%s ' % fname_outfilefile
cmd += 'ufafile=%s ' % fname_ufaevt
cmd += 'clfile=%s ' % fname_clevt
cmd += 'gtifile=%s ' % fname_gtifile
cmd += '>& %s ' % fname_log_nicermergeclean
print(cmd);os.system(cmd)

cmd = 'gzip %s;\n' % fname_ufaevt
cmd += 'gzip %s;\n' % fname_clevt
print(cmd);os.system(cmd)

fname_ufaevt += '.gz'
fname_clevt += '.gz'

# ==============================
# 5. make_bgdspec_mit3C50.py
# ==============================
dir_speclc = '%s/speclc' % dir_proc
prefix = 'ni%s_%s' % (obsid,param['gtiexpr_basestr'])

cmd  = 'make_bgdspec_mit3C50.py '
cmd += '%s %s' % (fname_ufaevt,fname_clevt)
cmd += '--outdir %s ' % dir_speclc
cmd += '--prefix %s ' % prefix
cmd += '--tbin %.1f ' % param['mitbgd_tbin']
cmd += '--lctbin %.1f --lcemin %.1f --lcemax %.1f ' % (
        param['lc_tbin'],param['lc_emin'],param['lc_emax'])
print(cmd);os.system(cmd)









