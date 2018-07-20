#!/usr/bin/env python

__name__    = 'niprocess'
__author__  = 'Teru Enoto'
__version__ = '1.01'
__date__    = '2018 July 12'
__keyword__ = 'niprocess_v%s' % __version__

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
import glob
import argparse
import astropy.io.fits as pyfits 

PATH_TO_FPARAM_FILE = '%s/scripts/process/setup_yaml' % os.getenv('NICER_SOFT_PATH')
DEFAULT_FPARAM_FILE = '%s/pipeline_setup_default.yaml' % PATH_TO_FPARAM_FILE

# ==============================
# Get input parameters 
# ==============================
help_message = """
(example) %s obsid_path.lst 
""" % sys.argv[0]
parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('indir', metavar='indir',type=str,        
	help='Input directory name. The directory should be a single NICER observation directory, which in turn contains xti/{events_uf,events_cl,hk,auxil} subdirectories.')
parser.add_argument('--fparam', metavar='fparam',type=str,     
	default=DEFAULT_FPARAM_FILE,   
	help='yaml file for input parameters.')
parser.add_argument('--outdir', metavar='outdir',type=str,     
	default="niout",help='main output directory.')
parser.add_argument("-r", "--recreate",action="store_true",
	dest="flag_recreate",default=False,help='recreate flag.')
args = parser.parse_args()
print(args)


indir_path_list = glob.glob('%s/*/%s' % (os.getenv('NICER_PUBLIC_DATA_PATH'),args.indir))
if len(indir_path_list) != 1:
	sys.stderr.write('%s does not exist (or several candidates).\n' % args.indir)
	exit()
indir_path = indir_path_list[0]	
obsid = os.path.basename(args.indir)
print("indir_path: %s" % indir_path)
print("obsid: %s" % obsid)

if not os.path.exists(args.fparam):
	sys.stderr.write('error: file %s does not exist.\n' % args.fparam)
	quit()

param = yaml.load(open(args.fparam))
print(param)

# ==============================
# 0. Prepare directory 
# ==============================
obsiddir = '%s/%s' % (args.outdir,obsid)

if args.flag_recreate:
	cmd = 'rm -rf %s' % (obsiddir)
	print(cmd);os.system(cmd)

if not os.path.exists(obsiddir):
	cmd = 'mkdir -p %s/xti' % (obsiddir)
	print(cmd);os.system(cmd)

	pwd = os.getcwd()
	os.chdir(obsiddir)
	cmd  = 'cp -r %s/auxil .;' % indir_path
	cmd += 'ln -s %s/log .;' % indir_path
	print(cmd);os.system(cmd)
	os.chdir('xti')
	cmd  = 'ln -s %s/xti/hk .;' % indir_path
	cmd += 'ln -s %s/xti/event_uf .;' % indir_path
	cmd += 'mkdir event_cl;' 
	#cmd += 'cp -r %s/xti/event_cl .;' % indir_path
	print(cmd);os.system(cmd)
	os.chdir(pwd)

# ==============================
# 1. nicercal - apply standard NICER calibration
# ==============================
fname_mkffile = '%s/auxil/ni%s.mkf.gz' % (obsiddir,obsid)

fname_log_nicercal = '%s/xti/event_cl/1_nicercal.log' % obsiddir
fname_outfilefile = '%s/xti/event_cl/1_nicercal_ufalist.lis' % obsiddir
cmd = 'nicercal indir=%s outfilefile=%s >& %s' % (obsiddir,fname_outfilefile,fname_log_nicercal)
print(cmd);
if os.path.exists(fname_log_nicercal):
	print("...already processed. skipped.")
else:
	os.system(cmd)
	cmd = 'gzip %s' % fname_mkffile.replace('.gz','')
	print(cmd);os.system(cmd)

fname_ufamerge_evt = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (obsiddir,obsid)
cmd = 'nimpumerge infiles=@%s outfile=%s ' % (fname_outfilefile,fname_ufamerge_evt)
print(cmd);
if os.path.exists(fname_ufamerge_evt+'.gz'):
	print("...already processed. skipped.")
else:
	os.system(cmd)
	cmd = 'gzip %s' % fname_ufamerge_evt
	print(cmd);os.system(cmd)
fname_ufamerge_evt += '.gz'

# ==============================
# 1. niprefilter2 - create augmented NICER-specific filter file
# ==============================
fname_log_niprefilter2 = '%s/auxil/2_niprefilter2.log' % obsiddir

if os.path.exists(fname_log_niprefilter2):
	print("...niprefilter2. already processed. skipped.")
else:
	cmd = 'gunzip %s' % fname_mkffile
	print(cmd);os.system(cmd)

	cmd = 'niprefilter2 clobber=yes indir=%s infile=%s outfile=INFILE >& %s' % (obsiddir,fname_mkffile.replace('.gz',''),fname_log_niprefilter2)
	print(cmd);os.system(cmd)

	cmd = 'gzip %s' % fname_mkffile.replace('.gz','')
	print(cmd);os.system(cmd)

# ==============================
# 3. nimaketime - create standard screening good time intervals
# ==============================
workdir  = '%s/work/%s' % (obsiddir,param['gtiexpr_basestr'])
dir_proc = '%s/process' % workdir
dir_log  = '%s/log' % dir_proc
cmd = 'mkdir -p %s' % dir_log
print(cmd);os.system(cmd)

cmd = 'cp %s %s' % (args.fparam,dir_log)
print(cmd);os.system(cmd)

fname_log_nimaketime = '%s/3_nimaketime.log' % dir_log
fname_expr_nimaketime = '%s/3_nimaketime.expr' % dir_log
fname_gtifile = '%s/ni%s_%s.gti' % (dir_proc,obsid,param['gtiexpr_basestr'])
cmd = 'nimaketime infile=%s outfile=%s expr="%s" outexprfile="%s" >& %s' % (fname_mkffile, fname_gtifile, param['nimaketime_expr'], fname_expr_nimaketime, fname_log_nimaketime)
print(cmd);os.system(cmd)

# ==============================
# 4. nicermergeclean - combine per-MPU data and filter/screen
# ==============================
fname_log_nicermergeclean = '%s/4_nicermergeclean.log' % dir_log
fname_ufaevt = '%s/ni%s_0mpu7_ufa_%s.evt' % (dir_proc,obsid,param['gtiexpr_basestr'])
fname_clevt  = '%s/ni%s_0mpu7_cl_%s.evt' % (dir_proc,obsid,param['gtiexpr_basestr'])
#fname_ufaevt = '%s/xti/event_cl/ni%s_0mpu7_ufa_%s.evt' % (workdir,obsid,param['gtiexpr_basestr'])
#fname_clevt  = '%s/xti/event_cl/ni%s_0mpu7_cl_%s.evt' % (workdir,obsid,param['gtiexpr_basestr'])
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
cmd += '%s %s ' % (fname_ufaevt,fname_clevt)
cmd += '--outdir %s ' % dir_speclc
cmd += '--prefix %s ' % prefix
cmd += '--exclude %s ' % str(param['exclude_detid']).replace('[','').replace(']','').replace(' ','')
cmd += '--tbin %.1f ' % param['mitbgd_tbin']
cmd += '--lctbin %.1f --lcemin %.1f --lcemax %.1f ' % (
        param['mitbgd_ql_lc_tbin'],param['mitbgd_ql_lc_emin'],param['mitbgd_ql_lc_emax'])
print(cmd);os.system(cmd)

# ==============================
# 6. Plot energy-selected curve
# ==============================
if len(glob.glob('%s/product' % dir_speclc)) == 0:
	print("no gti. exit.")
	exit()
product_cl_evt = glob.glob('%s/product/ni%s_*_clscr_gtisel.evt' % (dir_speclc,obsid))[0]

flc_eband_list = []
for eband in param['lc_eband_list']:
	fname_flc_eband = '%s_%sto%skeV.flc' % (product_cl_evt.replace('.evt',''),str(eband[0]).replace('.','p'),str(eband[1]).replace('.','p'))
	cmd  = 'fselect_filter_energy.py '
	cmd += '-i %s ' % product_cl_evt
	cmd += '-o %s ' % fname_flc_eband
	cmd += '--emin %.2f ' % eband[0]
	cmd += '--emax %.2f ' % eband[1]
	print(cmd);os.system(cmd)
	flc_eband_list.append(fname_flc_eband)

product_bandflc = product_cl_evt.replace('.evt','_eband.flc')
product_bandflc_ps = product_bandflc.replace('.flc','_flc.ps')
cmd = 'lcurve %d <<EOF\n' % len(flc_eband_list)
for flc_name in flc_eband_list:
	cmd += '%s\n' % flc_name
cmd += '-\n'	
cmd += '%d\n' % param['mitbgd_ql_lc_tbin']
cmd += '%d\n' % param['lc_nbint']
cmd += '%s\n' % product_bandflc
cmd += 'yes\n'
cmd += '/xw\n'
cmd += '%d\n' % len(flc_eband_list)
i = 2
for eband in param['lc_eband_list']:
	cmd += 'lab y%d %.1f-%.1f keV\n' % (i,eband[0],eband[1])
	cmd += 'col %d on %d\n' % (i,i)
	cmd += 'lab rotate\n'
	i += 1 
cmd += 'la f %s (%s)\n' % (obsid,param['gtiexpr_basestr'])
cmd += 'lwid 5 on 1..100\n'
cmd += 'lab rotate\n'
cmd += 'hard %s/cps\n'	% product_bandflc_ps
cmd += 'exit'
print(cmd);os.system(cmd)
cmd = 'ps2pdf %s' % product_bandflc_ps
print(cmd);os.system(cmd)

cmd = 'rm -f %s; mv %s %s' % (product_bandflc_ps,
	os.path.basename(product_bandflc_ps).replace('.ps','.pdf'),
	os.path.dirname(product_bandflc_ps))
print(cmd);os.system(cmd)

# ==============================
# 7. Calculate Counts 
# ==============================
hdu = pyfits.open(product_cl_evt)
target  = hdu[0].header['OBJECT']
dateobs = hdu[0].header['DATE-OBS']
exposure = hdu[1].header['EXPOSURE']
title = '%s %s %s (%.1f s)' % (target, obsid,dateobs,exposure)

fname_read_xcm = glob.glob('%s/product/ni%s_*_read.xcm' % (dir_speclc,obsid))[0]
fname_residual_ps = product_cl_evt.replace('.evt','_spec_residual.ps')
cmd  = 'xspec<<EOF\n'
cmd += 'setplot device tmp/null\n'			
cmd += '@%s\n' % fname_read_xcm
cmd += 'setplot energy\n'
cmd += 'notice **-**\n'
cmd += 'ignore **-0.2 15.0-**\n'
cmd += 'setplot rebin 3 50 1 \n'
cmd += 'setplot rebin 3 50 2 \n'
cmd += 'ipl d\n'
cmd += 'r y -3.0 3.0\n'
cmd += 'r x 0.2 15.0\n'
cmd += 'lab 1 pos 0.2 0 loc 0 1 " " ls 3 col 4 \n'
cmd += 'lwid 5 \n'
cmd += 'lwid 5 on 1..100\n'
cmd += 'time off\n'
cmd += 'la t %s\n' % title 
cmd += 'la f %s (obs=black, green=bgd, red=residual)\n' % obsid 
cmd += 'hard %s/cps\n' % fname_residual_ps
cmd += 'exit\n'
cmd += 'data 2 none\n'
tmp_logfile_list = []
for eband in param['xspec_eband_list']:
	logfile = 'tmp_rate_%sto%skeV.log' % (str(eband[0]).replace('.','p'),str(eband[1]).replace('.','p'))
	cmd += 'notice **-**\n'
	cmd += 'log %s\n' % logfile
	cmd += 'ignore **-%.2f %.2f-**\n' % (eband[0],eband[1])
	cmd += 'show rate\n'
	cmd += 'log none\n'
	tmp_logfile_list.append(logfile)
cmd += 'exit\n'	
print(cmd);os.system(cmd)

cmd = 'ps2pdf %s' % fname_residual_ps
print(cmd);os.system(cmd)

cmd = 'rm -f %s; mv %s %s' % (fname_residual_ps,
	os.path.basename(fname_residual_ps).replace('.ps','.pdf'),
	os.path.dirname(fname_residual_ps))
print(cmd);os.system(cmd)

xspec_rate_werror_list = []
for logfile in tmp_logfile_list:
	i = tmp_logfile_list.index(logfile)
	emin = param['xspec_eband_list'][i][0]
	emax = param['xspec_eband_list'][i][1]
	for line in open(logfile):
		cols = line.split()
		if cols[0] == '#Net':
			rate = float(cols[6])
			rate_error = float(cols[8])
			print(emin,emax,rate,rate_error)
			xspec_rate_werror_list.append([rate,rate_error])
	cmd = 'rm -f %s\n' % logfile
	print(cmd);os.system(cmd)


param['xspec_rate_werror_list'] = xspec_rate_werror_list
fname_result_yaml = fname_read_xcm.replace('_read.xcm','_out.yaml')
with open(fname_result_yaml, 'w') as yaml_file:
    yaml.dump(param, yaml_file, default_flow_style=True)

exit()





