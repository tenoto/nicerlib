#!/usr/bin/env python

import os 
import sys 
import yaml 
import glob 
import argparse
import numpy as np 
import pandas as pd 
import astropy.io.fits as pyfits 

help_message = """
(example) %s --indir indir --fparam fparam
""" % sys.argv[0]
parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('indir', metavar='indir',type=str,
	help='Input directory.') 
parser.add_argument('obsid_base', metavar='obsid_base',type=str,
	help='ObsID') 
parser.add_argument('outdir', metavar='outdir',type=str,
	help='output directory.') 
parser.add_argument('fname_input_xcm', metavar='fname_input_xcm',type=str,
	help='output directory.') 
parser.add_argument('model_type', metavar='model_type',type=str,
	help='model type (bb,pl).') 
parser.add_argument('--fit_emin', metavar='fit_emin',type=float,        
	help='fitting energy min',
	default=0.4)
parser.add_argument('--fit_emax', metavar='fit_emax',type=float,        
	help='fitting energy max',
	default=8.0)
parser.add_argument('--subtitle', metavar='subtitle',type=str,        
	help='subtitle',
	default="") 
parser.add_argument('--min_significance', metavar='min_significance',
	type=float, help='binning minimum significance',
	default=5)
parser.add_argument('--max_bins', metavar='max_bins',
	type=int, help='binning max bins',
	default=120)
parser.add_argument('--plot_ymin', metavar='plot_ymin',type=float,        
	help='plot_ymin',default=3e-4)
parser.add_argument('--plot_ymax', metavar='plot_ymax',type=float,        
	help='plot_ymax',default=300)
parser.add_argument('--plot_xmin', metavar='plot_xmin',type=float,        
	help='plot_xmin',default=0.2)
parser.add_argument('--plot_xmax', metavar='plot_xmax',type=float,        
	help='plot_xmax',default=10.0)
args = parser.parse_args()
print(args)

indir = args.indir 
obsid_base = args.obsid_base
dir_main = args.outdir
fname_input_xcm = args.fname_input_xcm
model_type = args.model_type

flag_execute = True

NULL_VALUE = -999
RCHI2_THRESHOLD = 4.0

# ========================
# functions
# ========================
def bin_spec(src_pha,bgd_pha,
	min_significance=5,max_bins=120,
	emin=0.2,emax=15.0):
	basename = '%s_grp' % os.path.splitext(os.path.basename(src_pha))[0]
	grp_qdp = '%s.qdp' % basename
	grp_pha = '%s.pha' % basename
	cmd  = 'xspec <<EOF\n'
	cmd += 'data 1 %s\n' % src_pha
	cmd += 'back 1 %s\n' % bgd_pha
	cmd += 'setplot energy\n'
	cmd += 'ignore 1:**-%.2f %.2f-**\n' % (emin,emax)
	cmd += 'setplot rebin %d %d\n' % (min_significance,max_bins)
	cmd += 'setplot channel\n'
	cmd += 'iplot data\n'
	cmd += 'we %s\n' % basename
	cmd += 'exit\n'
	cmd += 'exit\n'
	print(cmd);
	if flag_execute:
		os.system(cmd)

	cmd  = 'fgrppha.py %s %s' % (src_pha,grp_qdp)
	print(cmd);
	if flag_execute:
		os.system(cmd)

	cmd  = 'mv %s %s;' % (grp_pha,os.path.dirname(src_pha))
	cmd += 'rm -f %s.{qdp,pco};' % basename
	print(cmd);
	if flag_execute:
		os.system(cmd)

	grp_pha = '%s/%s' % (os.path.dirname(src_pha),grp_pha)
	return grp_pha

def fit_spec(src_grp_pha,bgd_pha,model_xcm,n_err_list=[1,2,5],subtitle=None,
	emin=0.2,emax=4.0,xmin=0.2,xmax=6,
	ymin=3e-3,ymax=300,y2min=-15,y2max=15,model_type='pl'):

	hdu = pyfits.open(src_grp_pha)
	date_obs   = hdu[1].header['DATE-OBS']
	obsid      = hdu[1].header['OBS_ID']
	objectname = hdu[1].header['OBJECT']
	exposure   = float(hdu[1].header['EXPOSURE'])
	title      = '%s %s %s (%.1f s)' % (objectname, obsid, date_obs, exposure)
	subtitle += 'fit at %.1f-%.1f keV (%s)' % (emin,emax,model_xcm)

	fname_fit_basename = '%s_fit' % (os.path.splitext(os.path.basename(src_grp_pha))[0])
	fname_fit_log = '%s/%s_fit.log' % (dir_log,os.path.splitext(os.path.basename(src_grp_pha))[0])
	fname_fit_xcm = fname_fit_log.replace('.log','.xcm')
	fname_fit_ps = fname_fit_log.replace('.log','.ps')
	cmd  = 'xspec <<EOF\n'
	cmd += 'data 1 %s\n' % src_grp_pha
	cmd += 'back 1 %s\n' % bgd_pha
	cmd += 'setplot energy\n'
	cmd += 'ignore 1:**-%0.2f %0.2f-**\n' % (emin,emax)
	cmd += '@%s\n' % model_xcm
	cmd += 'query yes\n'
	if model_type == 'bb':
		cmd += 'newpa 5 0 -1\n'
		cmd += 'freeze 4 5\n'	
		cmd += 'fit\n'
		cmd += 'freeze 1 2 3 \n'
		cmd += 'thaw 4 5\n'
		cmd += 'fit\n'
		cmd += 'thaw 2 3 \n'
		cmd += 'fit\n'
		cmd += 'thaw 1 \n' 
		cmd += 'fit\n'
	elif model_type == 'pl':
		cmd += 'fit\n'
		#cmd += 'setplot add\n'
	cmd += 'log %s\n' % fname_fit_log 
	cmd += 'show rate\n'
	cmd += 'show fit\n'
	cmd += 'show pa\n'
	cmd += 'freeze 1\n'
	cmd += 'fit\n'
	cmd += 'flux 1.0 10.0 err 100 68.3\n'
	cmd += 'flux 2.0 10.0 err 100 68.3\n'
	cmd += 'flux 0.2 1.0 err 100 68.3\n'	
	cmd += 'flux 0.4 6.0 err 100 68.3\n'		
	cmd += 'thaw 1\n'
	cmd += 'fit \n'
	for n_err in n_err_list:
		cmd += 'err 1.0 %d\n' % n_err
	cmd += 'log none\n'
	cmd += 'save all %s\n' % fname_fit_xcm
	cmd += 'iplot ld del\n'
	cmd += 'time off\n'
	cmd += 'la t %s \n' % title 
	if subtitle != None:
		cmd += 'la f %s\n' % subtitle
	cmd += 'lwid 5 \n'
	cmd += 'lwid 5 on 1..100 \n'	
	cmd += 'csize 1.1\n'
	cmd += 'lab pos y 2.8\n'
	cmd += 'r x %.1f %.1f\n' % (xmin,xmax)
	cmd += 'r y %.1e %.1e\n' % (ymin,ymax)
	cmd += 'r y2 %.1f %.1f\n' % (y2min,y2max)
	cmd += 'col 2 on 2\n'
	cmd += 'win 2\n'
	cmd += 'LAB  2 COL 2 LIN 0 100 JUS Lef POS 0.200000003 0 " "\n'
	cmd += 'hard %s/cps\n' % fname_fit_ps 
	cmd += 'we %s\n' % fname_fit_basename
	cmd += 'exit\n'
	cmd += 'exit\n'	
	print(cmd);
	if flag_execute:
		os.system(cmd)

	cmd = 'ps2pdf %s' % fname_fit_ps
	print(cmd);
	if flag_execute:
		os.system(cmd)
	fname_fit_pdf = os.path.basename(fname_fit_ps).replace('.ps','.pdf')
	cmd = 'mv %s %s; rm -f %s' % (fname_fit_pdf,os.path.dirname(fname_fit_ps),fname_fit_ps)
	print(cmd);
	if flag_execute:
		os.system(cmd)

	"""
	cmd  = 'echo "NO NO NO NO NO" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.35 0.15 3.47e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.75 0.25 1.40e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "1.50 0.50 6.85e-02 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "3.50 1.50 5.93e-02 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "7.50 2.50 3.90e-02 0" >> %s.qdp;\n' % fname_fit_basename 
	cmd += 'echo "12.50 2.50 2.51e-02 0" >> %s.qdp;\n' % fname_fit_basename
	print(cmd);os.system(cmd)
	"""

	# 5-sigma and 3-sigma uncertainty 
	cmd  = 'echo "NO NO NO NO NO" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.35 0.15 7.51e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.75 0.25 4.76e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "1.50 0.50 4.20e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "3.50 1.50 1.49e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "7.50 2.50 8.80e-02 0" >> %s.qdp;\n' % fname_fit_basename 
	cmd += 'echo "12.50 2.50 4.86e-02 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "NO NO NO NO NO" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.35 0.15 4.17e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "0.75 0.25 2.81e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "1.50 0.50 2.53e-01 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "3.50 1.50 8.27e-02 0" >> %s.qdp;\n' % fname_fit_basename
	cmd += 'echo "7.50 2.50 5.01e-02 0" >> %s.qdp;\n' % fname_fit_basename 
	cmd += 'echo "12.50 2.50 2.83e-02 0" >> %s.qdp;\n' % fname_fit_basename	
	print(cmd);os.system(cmd)

	cmd  = 'qdp %s.qdp<<EOF\n' % fname_fit_basename
	cmd += '/xw\n'
	cmd += '@%s.pco\n' % fname_fit_basename
	cmd += 'line step on 6\n' 
	cmd += 'col 8 on 6\n'
	cmd += 'line step on 5\n' 
	cmd += 'col 9 on 5\n'	
	cmd += 'hard %s_lim.ps/cps\n' % fname_fit_basename
	cmd += 'exit\n'
	print(cmd);os.system(cmd)

	cmd = 'ps2pdf %s_lim.ps' % fname_fit_basename
	print(cmd);os.system(cmd)	
	cmd = 'rm -f %s_lim.ps' % fname_fit_basename
	print(cmd);os.system(cmd)	

	cmd = 'mv %s.{qdp,pco} %s_lim.pdf %s' % (fname_fit_basename,fname_fit_basename,os.path.dirname(fname_fit_ps))
	print(cmd);os.system(cmd)

	return fname_fit_log

def get_flux(logfile,emin,emax):
	flag_found = False
	for line in open(logfile):
		cols = line.split() 
		if flag_found:
			flux_min = float(cols[6].replace('(',''))
			flux_max = float(cols[8].replace(')',''))
			break 
		if not "Model Flux" in line:
			continue
		emin_data = float(cols[8].replace('(',''))
		emax_data = float(cols[10])
		if emin == emin_data and emax == emax_data:
			print(line)
			flux = float(cols[5].replace('(',''))
			flag_found = True 
	flux_err_min = flux_min - flux
	flux_err_max = flux_max - flux 
	if not flag_found:
		print("no corresponding flux...")
	return flux, flux_err_min, flux_err_max

def get_rate(logfile):
	for line in open(logfile):
		cols = line.split()
		if "#Net count rate" in line:
			rate = float(cols[6])
			rate_err = float(cols[8])
			break
	return rate, rate_err

def get_chisqure(logfile):
	for line in open(logfile):
		cols = line.split()
		if "#Test statistic : Chi-Squared =" in line:
			chi2 = float(cols[5])
		if "# Reduced chi-squared =" in line:
			print(line)
			rchi2 = float(cols[4])
			dof = int(cols[6])
		if "# Null hypothesis probability =" in line:
			print(line)
			prob = float(cols[5])
	return rchi2, chi2, dof, prob

def get_parameter(logfile,npar):
	flag_show_pa = False
	for line in open(logfile):
		cols = line.split()
		if flag_show_pa and len(cols) > 1:
			if cols[1] == str(npar):
				value = float(cols[-3])
				break
		if "#Parameters defined:" in line:
			flag_show_pa = True

	rchi2, chi2, dof, prob = get_chisqure(logfile)
	if rchi2 > 2.0:
		err_min = NULL_VALUE
		err_max = NULL_VALUE
	else:
		flag_err = False
		for line in open(logfile):
			cols = line.split()
			if flag_err and len(cols) == 5:
				if cols[1] == str(npar) and cols[4][0] == '(' and cols[4][-1] == ')':
					print(line)
					tmp1,tmp2 = cols[4].split(',')
					err_min = float(tmp1.replace('(',''))
					err_max = float(tmp2.replace(')',''))
					break
			if "# Parameter   Confidence Range (1)" in line:
				flag_err = True

	return value, err_min, err_max

def get_MJD_UTC(date_utc):
	cmd  = 'rm -f tmp.log;'
	cmd += 'nitimeconv.py %s ' % date_utc
	cmd += '-f isot -s utc > tmp.log'
	print(cmd);os.system(cmd)

	for line in open('tmp.log'):
		cols = line.split()
		if len(cols) <= 1:
			continue 
		if cols[0] == 'MJD_UTC':
			mjd_utc = float(cols[2])
			break 
	cmd = 'rm -f tmp.log;'
	return mjd_utc

# ========================
# copy files 
# ========================
dir_input = '%s/data' % dir_main
dir_pha = '%s/pha' % dir_input
dir_pdf = '%s/pdf' % dir_input
dir_log = '%s/fit' % dir_input
cmd = 'rm -rf %s; mkdir -p %s %s %s' % (dir_main,dir_pha,dir_pdf,dir_log)
print(cmd);
if flag_execute:
	os.system(cmd)

for obsid_path in glob.glob('%s/%s*' % (indir,obsid_base)):
	for pha in glob.glob('%s/proc/speclc/prod/*.pha' % obsid_path):
		cmd = 'cp %s %s' % (pha,dir_pha)
		print(cmd);
		if flag_execute:
			os.system(cmd)
	for pdf in glob.glob('%s/proc/speclc/prod/*.pdf' % obsid_path):
		cmd = 'cp %s %s' % (pdf,dir_pdf)
		print(cmd);
		if flag_execute:
			os.system(cmd)

# ========================
# make fit 
# ========================
i = 0 
row_list = []
for src_pha in glob.glob('%s/*gtisel.pha' % dir_pha):
	hdu = pyfits.open(src_pha)

	date_obs = hdu[1].header['DATE-OBS']
	date_end = hdu[1].header['DATE-END']	
	obsid    = hdu[1].header['OBS_ID']
	exposure = float(hdu[1].header['EXPOSURE'])
	bgd_pha  = glob.glob('%s/ni%s_*_BGMod_3C50.pha' % (dir_pha,obsid))[0]

	mjd_utc_start = get_MJD_UTC(date_obs)
	mjd_utc_stop = get_MJD_UTC(date_end)	
	mjd_utc_center = 0.5 * (mjd_utc_start+mjd_utc_stop)
	mjd_utc_width = 0.5 * (-mjd_utc_start+mjd_utc_stop)	

	src_grp_pha = bin_spec(src_pha,bgd_pha,
		min_significance=args.min_significance,max_bins=args.max_bins)
	fname_fit_log = fit_spec(src_grp_pha,bgd_pha,fname_input_xcm,
		emin=args.fit_emin, emax=args.fit_emax, 
		subtitle=args.subtitle,ymin=args.plot_ymin,ymax=args.plot_ymax,
		xmin=args.plot_xmin,xmax=args.plot_xmax,model_type=model_type)

	rate, rate_err = get_rate(fname_fit_log)
	rchi2, chi2, dof, prov = get_chisqure(fname_fit_log)

	if rchi2 < RCHI2_THRESHOLD: 
		f1to10, f1to10_err_min, f1to10_err_max = get_flux(fname_fit_log,1.0,10.0)
		f2to10, f2to10_err_min, f2to10_err_max = get_flux(fname_fit_log,2.0,10.0)	
		f0p2to1, f0p2to1_err_min, f0p2to1_err_max = get_flux(fname_fit_log,0.2,1.0)
		f0p4to6, f0p4to6_err_min, f0p4to6_err_max = get_flux(fname_fit_log,0.4,6.0)	
		if model_type == 'bb':
			nh, nh_err_min, nh_err_max = get_parameter(fname_fit_log,1)
			kT1, kT1_err_min, kT1_err_max = get_parameter(fname_fit_log,2)
			norm1, norm1_err_min, norm1_err_max = get_parameter(fname_fit_log,3)
		if model_type == 'pl':
			nh, nh_err_min, nh_err_max = get_parameter(fname_fit_log,1)
			gamma, gamma_err_min, gamma_err_max = get_parameter(fname_fit_log,2)
			norm1, norm1_err_min, norm1_err_max = get_parameter(fname_fit_log,5)			
	else:
		f1to10, f1to10_err_min, f1to10_err_max = np.nan, np.nan, np.nan
		f2to10, f2to10_err_min, f2to10_err_max = np.nan, np.nan, np.nan
		f0p2to1, f0p2to1_err_min, f0p2to1_err_max = np.nan, np.nan, np.nan
		f0p4to6, f0p4to6_err_min, f0p4to6_err_max = np.nan, np.nan, np.nan
		if model_type == 'bb':
			nh, nh_err_min, nh_err_max = np.nan, np.nan, np.nan
			kT1, kT1_err_min, kT1_err_max = np.nan, np.nan, np.nan
			norm1, norm1_err_min, norm1_err_max = np.nan, np.nan, np.nan
		if model_type == 'pl':
			nh, nh_err_min, nh_err_max = np.nan, np.nan, np.nan
			gamma, gamma_err_min, gamma_err_max = np.nan, np.nan, np.nan
			norm1, norm1_err_min, norm1_err_max = np.nan, np.nan, np.nan

	if model_type == 'pl':
		row_list.append([
			date_obs,
			mjd_utc_start,
			mjd_utc_stop,
			mjd_utc_center,	
			mjd_utc_width,				
			obsid,
			exposure,
			rate,
			rate_err,
			args.fit_emin,
			args.fit_emax,
			rchi2,
			chi2,
			dof,
			prov,
			f1to10, f1to10_err_min, f1to10_err_max,
			f2to10, f2to10_err_min, f2to10_err_max,
			f0p2to1, f0p2to1_err_min, f0p2to1_err_max,
			f0p4to6, f0p4to6_err_min, f0p4to6_err_max,
			nh, nh_err_min, nh_err_max,
			gamma, gamma_err_min, gamma_err_max,
			norm1, norm1_err_min, norm1_err_max, 
			src_pha,
			bgd_pha
			])
	else:
		row_list.append([
			date_obs,
			mjd_utc_start, 
			mjd_utc_stop,	
			mjd_utc_center,	
			mjd_utc_width,			
			obsid,
			exposure,
			rate,
			rate_err,
			args.fit_emin,
			args.fit_emax,
			rchi2,
			chi2,
			dof,
			prov,
			f1to10, f1to10_err_min, f1to10_err_max,
			f2to10, f2to10_err_min, f2to10_err_max,
			f0p2to1, f0p2to1_err_min, f0p2to1_err_max,
			f0p4to6, f0p4to6_err_min, f0p4to6_err_max,
			src_pha,
			bgd_pha
			])		

	i += 1 
	#if i > 3:
	#	break

print(row_list)

fname_main_log = '%s/%s_fit.csv' % (dir_main,args.outdir)
df = pd.DataFrame(row_list,
	columns=[
	'DATE-OBS',
	'MJD_UTC_START',
	'MJD_UTC_STOP',	
	'MJD_UTC_CENTER',		
	'MJD_UTC_WIDTH',			
	'OBS_ID',
	'EXPOSURE(s)',
	'rate',
	'rate_err',
	'fit_emin',
	'fit_emax',
	'rchi2',
	'chi2',
	'dof',
	'prov',
	'f1to10', 'f1to10_err_min', 'f1to10_err_max',
	'f2to10', 'f2to10_err_min', 'f2to10_err_max',
	'f0p2to1', 'f0p2to1_err_min', 'f0p2to1_err_max',
	'f0p4to6', 'f0p4to6_err_min', 'f0p4to6_err_max',
	'nh', 'nh_err_min', 'nh_err_max',
	'gamma', 'gamma_err_min', 'gamma_err_max',
	'norm1', 'norm1_err_min', 'norm1_err_max', 	
	'src_pha',
	'bgd_pha'
	])		
df.to_csv(fname_main_log)

"""
# ========================
# make fit 
# ========================
fname_main_log = '%s/README_%s.txt' % (dir_main,args.outdir)
f_main_log = open(fname_main_log,'w')
dump  = '# %s' % args.outdir
dump += '# ObsID DATE-OBS Exp(s)'
dump += 'rate, rate_err  '
dump += 'rchi2, chi2, dof, prov  '
dump += 'f0p4to6,f0p4to6_err_min,f0p4to6_err_max  '
dump += 'nh, nh_err_min, nh_err_max  '
dump += 'kT1, kT1_err_min, kT1_err_max  '
dump += 'norm1, norm1_err_min, norm1_err_max  '
dump += 'f1to10, f1to10_err_min, f1to10_err_max  '
dump += 'f2to10,f2to10_emin,f2to10_err_max  '
dump += 'f0p2to1, f0p2to1_err_min, f0p2to1_err_max  '
dump += '\n'
f_main_log.write(dump)
for src_pha in glob.glob('%s/*gtisel.pha' % dir_pha):
	print(src_pha)
	hdu = pyfits.open(src_pha)
	date_obs = hdu[1].header['DATE-OBS']
	obsid = hdu[1].header['OBS_ID']
	exposure = float(hdu[1].header['EXPOSURE'])
	bgd_pha = glob.glob('%s/ni%s_*_BGMod_3C50.pha' % (dir_pha,obsid))[0]
	print(src_pha, bgd_pha)
	src_grp_pha = bin_spec(src_pha,bgd_pha,
		min_significance=args.min_significance,max_bins=args.max_bins)
	fname_fit_log = fit_spec(src_grp_pha,bgd_pha,fname_input_xcm,
		emin=args.fit_emin, emax=args.fit_emax, 
		subtitle=args.subtitle,ymin=args.plot_ymin,ymax=args.plot_ymax,
		xmin=args.plot_xmin,xmax=args.plot_xmax)

	rate, rate_err = get_rate(fname_fit_log)
	rchi2, chi2, dof, prov = get_chisqure(fname_fit_log)
	dump_fit = '%.3f %.3f  ' % (rate, rate_err)
	dump_fit += '%.2f %.2f %d %.2e  ' % (rchi2, chi2, dof, prov)		

	if rchi2 < RCHI2_THRESHOLD: 
		f1to10, f1to10_err_min, f1to10_err_max = get_flux(fname_fit_log,1.0,10.0)
		f2to10, f2to10_err_min, f2to10_err_max = get_flux(fname_fit_log,2.0,10.0)	
		f0p2to1, f0p2to1_err_min, f0p2to1_err_max = get_flux(fname_fit_log,0.2,1.0)
		f0p4to6, f0p4to6_err_min, f0p4to6_err_max = get_flux(fname_fit_log,0.4,6.0)	
		nh, nh_err_min, nh_err_max = get_parameter(fname_fit_log,1)
		kT1, kT1_err_min, kT1_err_max = get_parameter(fname_fit_log,2)
		norm1, norm1_err_min, norm1_err_max = get_parameter(fname_fit_log,3)
		dump_fit += '%.3e %.3e %.3e  ' % (f0p4to6,f0p4to6_err_min,f0p4to6_err_max)
		dump_fit += '%.3f %.3f %.3f  ' % (nh, nh_err_min, nh_err_max)
		dump_fit += '%.3f %.3f %.3f  ' % (kT1, kT1_err_min, kT1_err_max)
		dump_fit += '%.3f %.3f %.3f  ' % (norm1, norm1_err_min, norm1_err_max)
		dump_fit += '%.3e %.3e %.3e ' % (f1to10, f1to10_err_min, f1to10_err_max)
		dump_fit += '%.3e %.3e %.3e ' % (f2to10,f2to10_err_min,f2to10_err_max)	
		dump_fit += '%.3e %.3e %.3e ' % (f0p2to1, f0p2to1_err_min, f0p2to1_err_max)	
	else:
		dump_fit += '- - -   ' 
		dump_fit += '- - -   ' 
		dump_fit += '- - -   ' 
		dump_fit += '- - -   ' 
		dump_fit += '- - -   ' 
		dump_fit += '- - -   '
		dump_fit += '- - -   ' 

	dump = '%s  %s  %.1f  ' % (obsid,date_obs,exposure)
	dump += dump_fit
	dump += '\n'
	f_main_log.write(dump)

	fname_fit_yaml = fname_fit_log.replace('.log','.yaml')
	dict_param = {}
	dict_param['DATE-OBS'] = date_obs
	dict_param['OBSID'] = obsid
	dict_param['EXPOSURE'] = exposure 
	dict_param['src_pha'] = src_pha
	dict_param['src_grp_pha'] = src_grp_pha
	dict_param['bgd_pha'] = bgd_pha
	dict_param['rate'] = rate
	dict_param['rate_err'] = rate_err
	dict_param['rchi2'] = rchi2
	dict_param['chi2'] = chi2
	dict_param['dof'] = dof
	dict_param['prov'] = prov

	dict_param['f1to10'] = f1to10
	dict_param['f1to10_err_min'] = f1to10_err_min
	dict_param['f1to10_err_max'] = f1to10_err_max
	dict_param['f2to10'] = f2to10
	dict_param['f2to10_err_min'] = f2to10_err_min
	dict_param['f2to10_err_max'] = f2to10_err_max	
	dict_param['f0p2to1'] = f0p2to1
	dict_param['f0p2to1_err_min'] = f0p2to1_err_min
	dict_param['f0p2to1_err_max'] = f0p2to1_err_max	
	dict_param['f0p4to6'] = f0p4to6
	dict_param['f0p4to6_err_min'] = f0p4to6_err_min
	dict_param['f0p4to6_err_max'] = f0p4to6_err_max		
	dict_param['nh'] = nh
	dict_param['nh_err_min'] = nh_err_min
	dict_param['nh_err_max'] = nh_err_max	
	dict_param['kT1'] = kT1
	dict_param['kT1_err_min'] = kT1_err_min
	dict_param['kT1_err_max'] = kT1_err_max		
	dict_param['norm1'] = norm1
	dict_param['norm1_err_min'] = norm1_err_min
	dict_param['norm1_err_max'] = norm1_err_max

	f_fit_yaml = open(fname_fit_yaml,'w')
	f_fit_yaml.write(yaml.dump(dict_param,default_flow_style=False))
	f_fit_yaml.close()

	rchi2 = NULL_VALUE
	chi2 = NULL_VALUE
	dof = NULL_VALUE
	prov = NULL_VALUE

	f1to10 = NULL_VALUE
	f1to10_err_min = NULL_VALUE
	f1to10_err_max = NULL_VALUE
	f2to10 = NULL_VALUE
	f2to10_err_min = NULL_VALUE
	f2to10_err_max = NULL_VALUE
	f0p2to1 = NULL_VALUE
	f0p2to1_err_min = NULL_VALUE
	f0p2to1_err_max = NULL_VALUE
	nh = NULL_VALUE
	nh_err_min = NULL_VALUE
	nh_err_max = NULL_VALUE
	kT1 = NULL_VALUE
	kT1_err_min = NULL_VALUE
	kT1_err_max = NULL_VALUE
	norm1 = NULL_VALUE
	norm1_err_min = NULL_VALUE
	norm1_err_max = NULL_VALUE

f_main_log.close()
"""














