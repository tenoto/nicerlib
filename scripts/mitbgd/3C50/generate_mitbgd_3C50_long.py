#!/usr/bin/env python

__name__    = 'generate_mitbgd_3C50_long'
__author__  = 'Teru Enoto (for Rons model with Mikes application)'
__version__ = '1.00'
__date__    = '2018 June 6'

"""
see the document 'NICER XTI Background Model 3C50'
Description and User Instructions
(BGMod_3C50.tar)
v.3 2018 0531

with Mikes application 
"""

BGD_LIBRARY_PATH = "/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd/BGMod_3C50"
BG_GROUP_3C50_TABLE = "%s/bg_group_3C50.table" % BGD_LIBRARY_PATH
DAY_NOISE_3C50_TABLE = "%s/day_noise_3C50.table" % BGD_LIBRARY_PATH
NUM_OF_MPU = 7
NUM_OF_FPM = 8

import os
import sys
import subprocess
import numpy as np 
import astropy.io.fits as pyfits
import argparse
from datetime import datetime 

def get_lc_rate(lcfits):
	hdu = pyfits.open(lcfits)
	return np.mean(hdu['RATE'].data['RATE'])

class ExposureMap_IBGvsHREJ():
	def __init__(self):
		self.matrix_ibg52_hrej52 = np.zeros((6,5))

	def add(self,IBG52,HREJ52,exposure,NFPM):
		if IBG52 < 0.10:
			if HREJ52<0.15: 
				i = 0; j = 0;
			elif 0.15<=HREJ52<0.40:
				i = 0; j = 1;				
			elif 0.40<=HREJ52<1.00:
				i = 0; j = 2;								
			elif 1.00<=HREJ52<2.00:
				i = 0; j = 3;												
			elif 2.00<=HREJ52<5.00:
				i = 0; j = 4;																
		elif 0.10 <= IBG52 < 0.20:
			if HREJ52<0.18: 
				i = 1; j = 0;
			elif 0.18<=HREJ52<0.40:
				i = 1; j = 1;
			elif 0.40<=HREJ52<1.00:
				i = 1; j = 2;
			elif 1.00<=HREJ52<2.00:
				i = 1; j = 3; 
			elif 2.00<=HREJ52<5.00:
				i = 1; j = 4;
		elif 0.20 <= IBG52 < 0.40:
			if HREJ52<0.20: 
				i = 2; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 2; j = 1; 
			elif 0.40<=HREJ52<1.00:
				i = 2; j = 2 ;
			elif 1.00<=HREJ52<2.00:
				i = 2; j = 3;
			elif 2.00<=HREJ52<5.00:
				i = 2; j = 4;
		elif 0.40 <= IBG52 < 1.00:
			if HREJ52<0.20: 
				i = 3; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 3; j = 1;
			elif 0.40<=HREJ52<0.90:
				i = 3; j = 2;
			elif 0.90<=HREJ52<2.00:
				i = 3; j = 3;
			elif 2.00<=HREJ52<5.00:
				i = 3; j = 4;
		elif 1.00 <= IBG52 < 3.00:
			if HREJ52<0.20: 
				i = 4; j = 0;
			elif 0.20<=HREJ52<0.40:
				i = 4; j = 1;
			elif 0.40<=HREJ52<1.00:
				i = 4; j = 2; 
			elif 1.00<=HREJ52<2.30:
				i = 4; j = 3;
			elif 2.30<=HREJ52<5.00:
				i = 4; j = 4;
		elif 3.00 <= IBG52 < 10.00:
			if HREJ52<0.20: 
				i = 5; j = 0;
			elif 0.20<=HREJ52<0.55:
				i = 5; j = 1;
			elif 0.55<=HREJ52<1.00:
				i = 5; j = 2;
			elif 1.00<=HREJ52<2.00:
				i = 5; j = 3;
			elif 2.00<=HREJ52<5.00:				
				i = 5; j = 4;

		target_phaname = 'bg_group_3C50_ngt_%d%d.pha' % (i+1,j+1)
		for line2 in open(BG_GROUP_3C50_TABLE):
			cols2 = line2.split()
			phaname = cols2[0]
			if phaname == target_phaname:
				IBG_REF = float(cols2[1])
				break 

		weight = exposure * IBG52 / IBG_REF * NFPM / 52.0
		self.matrix_ibg52_hrej52[i][j] += weight		
		return i,j

	def show(self):
		print(self.matrix_ibg52_hrej52)

	def get_grand_exposure(self):
		return sum(sum(self.matrix_ibg52_hrej52))

	def get_normalized_exposuremap(self):
		total_exposure = self.get_grand_exposure()
		return self.matrix_ibg52_hrej52 / total_exposure

class ExposureMap_NZ():
	def __init__(self):
		self.hist_nz52 = np.zeros(15)

	def add(self,NZ52,exposure,NFPM):
		if NZ52 < 190.0: 
			k = 0
		elif 190.0 <= NZ52 < 210.0:
			k = 1
		elif 210.0 <= NZ52 < 250.0:			
			k = 2
		elif 250.0 <= NZ52 < 300.0:
			k = 3
		elif 300.0 <= NZ52 < 400.0:
			k = 4
		elif 400.0 <= NZ52 < 500.0:
			k = 5
		elif 500.0 <= NZ52 < 600.0:
			k = 6
		elif 600.0 <= NZ52 < 700.0:
			k = 7
		elif 700.0 <= NZ52 < 800.0:
			k = 8
		elif 800.0 <= NZ52 < 950.0:
			k = 9
		elif 950.0 <= NZ52 < 1100.0:
			k = 10
		elif 1100.0 <= NZ52 < 1250.0:
			k = 11
		elif 1250.0 <= NZ52 < 1400.0:
			k = 12
		elif 1400.0 <= NZ52 < 1600.0:
			k = 13
		elif 1600.0 <= NZ52:
			k = 14

		print(k)
		target_phaname = 'day_noise_nz%02d.pha' % k
		if k == 0 or k == 14:
			weight = 0.0
		else:
			for line2 in open(DAY_NOISE_3C50_TABLE):
				cols2 = line2.split()
				phaname = cols2[0]
				print(cols2)
				if phaname == target_phaname:
					NZ52_REF = float(cols2[1])
					break 
			weight = exposure * NZ52 / NZ52_REF * NFPM / 52.0
		self.hist_nz52[k] += weight 
		return k

	def show(self):
		print(self.hist_nz52)

	def get_grand_exposure(self):
		return sum(self.hist_nz52)

	def get_normalized_exposuremap(self):
		total_exposure = self.get_grand_exposure()
		if total_exposure == 0.0:
			return self.hist_nz52
		else:
			return self.hist_nz52 / total_exposure


cmd = 'which fhelp'
resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
if len(resp.stdout.readlines()) == 0:
	sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
	exit()

help_message = """
(example) %s obsid_path.lst 
""" % sys.argv[0]

parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('obsid_path', metavar='obsid_path',type=str,        
	help='input obsid_path or @input_obsid_path.lst for multiple obsids.')
parser.add_argument('--outdir',dest='outdir',action='store',
	help='output directory', default='mitbgd') 
parser.add_argument('--prefix',dest='prefix',action='store',
	help='prefix of the output files', default='mitbgd') 
parser.add_argument('--recreate', action='store_true', 
	help='flag to recreate of the output directory.') 
parser.add_argument('--exclude',dest='exclude',action='store',
	help='excluding detid (e.g., --excludedetid 14,34', default='14,34') 
parser.add_argument('--tbin',dest='tbin',action='store',
	help='GTI time scale (default 100 sec)', default=100.0) 
parser.add_argument('--tlimit',dest='tlimit',action='store',
	help='GTI time scale lower threshold (default 60 sec)', default=60.0) 
#parser.add_argument('--checkoutlier', action='store_true', 
#	help='flag to check outlier.') 
args = parser.parse_args()
print(args)

obsid_path_list = []
if args.obsid_path[0] is '@':
	for line in open(args.obsid_path[1:]):
		obsid_path_list.append(line.strip().strip('/'))
else:
	obsid_path_list.append(args.obsid_path)
print("obsid_path_list: %s" % obsid_path_list)

outdir = args.outdir
if args.recreate:
	cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
elif not os.path.exists(outdir):
	cmd = 'mkdir -p %s' % outdir
else:
	sys.stderr.write('output directory has already exited: %s\n' % outdir)
	quit()
print(cmd);os.system(cmd)

prefix = args.prefix

exclude_detid_list = []
for i in args.exclude.split(','):
	exclude_detid_list.append('%02d' % int(i))
print('exclude_detid_list: %s\n' % exclude_detid_list)

num_of_fpm = 56 - 4 - len(exclude_detid_list) 

obsid_list  = []
ufaevt_list = []
clevt_list  = []
for obsid_path in obsid_path_list:
	obsid = obsid_path.split('/')[-1]
	obsid_list.append(obsid)
	ufaevt    = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (obsid_path,obsid)
	ufaevt_gz = '%s.gz' % ufaevt
	if os.path.exists(ufaevt_gz):
		ufaevt_list.append(ufaevt_gz)
	elif os.path.exists(ufaevt):
		ufaevt_list.append(ufaevt)
	else:
		sys.stderr.write('ufa file does not exist: %s or %s \n' % (ufaevt,ufaevt_gz))
		quit()

	clevt    = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (obsid_path,obsid)
	clevt_gz = '%s.gz' % clevt
	if os.path.exists(clevt_gz):
		clevt_list.append(clevt_gz)
	elif os.path.exists(clevt):
		clevt_list.append(clevt)
	else:
		sys.stderr.write('ufa file does not exist: %s or %s \n' % (clevt,clevt_gz))
		quit()		
print("ufaevt list: %s" % ufaevt_list)
print("clevt list: %s" % clevt_list)

### GTI 
outdir_gti = '%s/gti' % outdir
cmd = 'rm -rf %s; mkdir -p %s' % (outdir_gti,outdir_gti)
print(cmd);os.system(cmd)

gti_list = []
f_gti_list = '%s/%s_gti.lst' % (outdir_gti,prefix)
f = open(f_gti_list,'w')
for clevt in clevt_list:
	obsid = os.path.basename(clevt).split('_')[0].split('ni')[-1]
	hdu = pyfits.open(clevt)
	print(clevt,obsid)
	for row in hdu['GTI'].data:
		gti_exp = row['STOP'] - row['START']
		if gti_exp < args.tlimit:
			message = '%.1f %.1f %.1f (skip)' % (row['START'],row['STOP'],gti_exp)
			print(message)
			continue
		ngti_subdiv = int(gti_exp/float(args.tbin))
		message = '%.1f %.1f %.1f %d' % (row['START'],row['STOP'],gti_exp,ngti_subdiv)
		print(message)		
		for i in range(ngti_subdiv+1):
			print(row['START'],row['STOP'],i,args.tbin)
			tstart = row['START'] + float(i) * float(args.tbin)
			tstop  = row['START'] + float(i+1) * float(args.tbin)
			if tstop > row['STOP']:
				tstop = row['STOP']
			tgti = tstop - tstart
			if tgti > args.tlimit:
				nid = len(gti_list)
				gti_list.append([nid,tstart,tstop,tgti,obsid])
				message = '  %d %.1f %.1f %.1f %s' % (nid,tstart,tstop,tgti,obsid)
				f.write('%.1f %.1f %d %s %.1f\n' % (tstart,tstop,nid,obsid,tgti))
			else:
				message = '  -- %.1f %.1f %.1f (skip)' % (tstart,tstop,tgti)
			print(message)
f.close()			
print('GTI: %s' % gti_list)

f_ufaevt_list = '%s/%s_ufaevt.lst' % (outdir,prefix)
f = open(f_ufaevt_list,'w')
for ufaevt in ufaevt_list:
	f.write(ufaevt+'\n')
f.close()

f_clevt_list = '%s/%s_clevt.lst' % (outdir,prefix)
f = open(f_clevt_list,'w')
for clevt in clevt_list:
	f.write(clevt+'\n')
f.close()

# 4.1 Make trimmed, slow-chain ufa files to support down-stream queries 
f_ufaslow = '%s/%s_ufaslow.evt' % (outdir,prefix)
cmd  = 'niextract-events '
cmd += 'filename="@%s[PI=-10:1800,EVENT_FLAGS=bxxx1x000]" ' % f_ufaevt_list
cmd += 'eventsout="%s" gti=GTI ' % f_ufaslow
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_ufaslow
print(cmd);os.system(cmd)
f_ufaslow += '.gz'
print(f_ufaslow)

# 4.2 Make trumpet-selected ufa files to screen FPMs & evaluate IBG, NZ
f_ufatsel = '%s/%s_ufatsel.evt' % (outdir,prefix)
cmd  = 'fselect %s %s ' % (f_ufaslow,f_ufatsel)
cmd += '"ISNULL(PI_RATIO) || PI_RATIO < (1.1+120/PI)" '
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_ufatsel
print(cmd);os.system(cmd)
f_ufatsel += '.gz'
print(f_ufatsel)

# 4.3 Make hatchet-rejected ufa files to evaluate HREJ
f_ufahrej = '%s/%s_ufahrej.evt' % (outdir,prefix)
cmd  = 'fselect %s %s ' % (f_ufaslow,f_ufahrej)
cmd += '"PI_RATIO >= 1.54" '
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_ufahrej
print(cmd);os.system(cmd)
f_ufahrej += '.gz'
print(f_ufahrej)

"""
# 4.3 4.4 FPM Screening: Event Lists per FPM
if args.checkoutlier:

	outdir_detid = '%s/detid' % outdir
	cmd = 'rm -rf %s; mkdir -p %s' % (outdir_detid,outdir_detid)
	print(cmd);os.system(cmd)

	for mpu in range(0,NUM_OF_MPU):
		for fpm in range(0,NUM_OF_FPM):
			detid = "%02d" % (10*mpu+fpm)
			f_ufatsel_detid = '%s/%s_ufatsel_detid%s.evt' % (outdir_detid,prefix,detid)
			cmd  = 'fselect %s %s ' % (f_ufatsel,f_ufatsel_detid)
			cmd += '"DET_ID == %s" ' % detid
			print(cmd);os.system(cmd)

			cmd  = 'gzip %s' % f_ufatsel_detid
			print(cmd);os.system(cmd)
			f_ufatsel_detid += '.gz'
			print(f_ufatsel_detid)

	f_gti_list_detail = '%s/%s_gti_detail.lst' % (outdir_gti,prefix)
	flog = open(f_gti_list_detail,'w')
	for gti in gti_list:
		gtinum = gti[0]
		tstart = gti[1]
		tstop  = gti[2]
		texp   = gti[3]
		subdir = '%s/gti%05d' % (outdir_gti,gtinum)
		cmd = 'rm -rf %s; mkdir -p %s' % (subdir,subdir)
		print(cmd);os.system(cmd)

		f_gti_outlier = '%s/%s_gti%d.txt' % (subdir,prefix,gtinum)
		f = open(f_gti_outlier,'w')

		gti_nz_rate_list = []
		gti_inband_rate_list = []
		for mpu in range(0,NUM_OF_MPU):
			for fpm in range(0,NUM_OF_FPM):
				detid = "%02d" % (10*mpu+fpm)
				f_ufatsel_detid = '%s/%s_ufatsel_detid%s.evt' % (outdir_detid,prefix,detid)

				f_ufatsel_detid_gti_nz = '%s/%s_ufatsel_detid%s_nz.evt' % (subdir,prefix,detid)
				cmd  = 'niextract-events '
				cmd += 'filename="%s[PI=0:200,TIME=%.1f:%.1f]" ' % (f_ufatsel_detid,tstart,tstop)
				cmd += 'eventsout="%s" gti=GTI ' % f_ufatsel_detid_gti_nz
				print(cmd);os.system(cmd)

				hdu = pyfits.open(f_ufatsel_detid_gti_nz)
				nz_count = len(hdu['EVENTS'].data)
				nz_rate = float(nz_count) / texp 
				gti_nz_rate_list.append(nz_rate)

				f_ufatsel_detid_gti_inband = '%s/%s_ufatsel_detid%s_inband.evt' % (subdir,prefix,detid)
				cmd  = 'niextract-events '
				cmd += 'filename="%s[PI=40:1200,TIME=%.1f:%.1f]" ' % (f_ufatsel_detid,tstart,tstop)
				cmd += 'eventsout="%s" gti=GTI ' % f_ufatsel_detid_gti_inband
				print(cmd);os.system(cmd)

				hdu = pyfits.open(f_ufatsel_detid_gti_inband)
				inband_count = len(hdu['EVENTS'].data)
				inband_rate = float(inband_count) / texp 
				gti_inband_rate_list.append(inband_rate)

				dump = '%s %d %d %.4f %.4f' % (detid, mpu, fpm, nz_rate, inband_rate)
				f.write(dump+'\n')
				print(dump)
		print("gti_nz_rate_list: %s" % gti_nz_rate_list)
		print("gti_inband_rate_list: %s" % gti_inband_rate_list)
		f.close()	

		gti_nz_rate_list.remove(0.0)
		gti_inband_rate_list.remove(0.0)

		gti_nz_rate_list_sort = sorted(gti_nz_rate_list)
		gti_nz_rate_mean = np.mean(gti_nz_rate_list_sort[2:-2])
		gti_nz_rate_std = np.std(gti_nz_rate_list_sort[2:-2])	

		gti_inband_rate_list_sort = sorted(gti_inband_rate_list)
		gti_inband_rate_mean = np.mean(gti_inband_rate_list_sort[2:-2])
		gti_inband_rate_std = np.std(gti_inband_rate_list_sort[2:-2])		

		nz_outlier = []
		inband_outlier = []
		for line in open(f_gti_outlier):
			cols = line.split()
			detid = cols[0]
			mpu = int(cols[1])
			fpm = int(cols[2])
			nz_rate = float(cols[3])
			inband_rate = float(cols[4])

			if (nz_rate < gti_nz_rate_mean-4*gti_nz_rate_std) or (nz_rate > gti_nz_rate_mean+4*gti_nz_rate_std):
				nz_outlier.append(detid)
			if (inband_rate < gti_inband_rate_mean-4*gti_inband_rate_std) or (inband_rate > gti_inband_rate_mean+4*gti_inband_rate_std):
				inband_outlier.append(detid)

		dump  = '%d %.1f %.1f %.3f  ' % (gtinum, tstart, tstop, texp)
		dump += '%.4f %.4f %s  ' % (gti_nz_rate_mean, gti_nz_rate_std, nz_outlier)
		dump +=' %.4f %.4f %s  ' % (gti_inband_rate_mean, gti_inband_rate_std, inband_outlier) 
		flog.write(dump+'\n')
	flog.close()	
"""

# 4.5 FPM Screening the cl event lists prior to extracting spectra per GTI
f_clevt_merge = '%s/%s_clmerge.evt' % (outdir,prefix)
cmd  = 'ftmerge @%s %s' % (f_clevt_list,f_clevt_merge)
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_clevt_merge
print(cmd);os.system(cmd)
f_clevt_merge += '.gz'

# screen of detid 
f_clevt_screen = '%s/%s_clmerge_screen.evt' % (outdir,prefix)
if len(exclude_detid_list) == 0:
	cmd = 'cp %s %s ' % (f_clevt_merge, f_clevt_screen)
else:
	cmd  = 'fselect %s %s ' % (f_clevt_merge,f_clevt_screen)
	cmd += '"'
	for detid in exclude_detid_list:
		cmd += '(DET_ID != %s)' % (detid)
		if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
			cmd += ' && '
	cmd += '"'
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_clevt_screen
print(cmd);os.system(cmd)
f_clevt_screen += '.gz'

# 4.6 FPM Screening ufa event lists prior to extracting IBG, HREJ, NZ
# do 4.2 fselect command with addition "(DET_ID != 14) && (DET_ID != 34)"
# do 4.3 fselect command with addition "(DET_ID != 14) && (DET_ID != 34)"

f_ufatsel_screen = '%s/%s_ufatsel_screen.evt' % (outdir,prefix)
if len(exclude_detid_list) == 0:
	cmd = 'cp %s %s ' % (f_ufatsel, f_ufatsel_screen)
else:
	cmd  = 'fselect %s %s ' % (f_ufatsel,f_ufatsel_screen)
	cmd += '"'
	for detid in exclude_detid_list:
		cmd += '(DET_ID != %s)' % (detid)
		if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
			cmd += ' && '
	cmd += '"'
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_ufatsel_screen
print(cmd);os.system(cmd)
f_ufatsel_screen += '.gz'

f_ufahrej_screen = '%s/%s_ufahrej_screen.evt' % (outdir,prefix)
if len(exclude_detid_list) == 0:
	cmd = 'cp %s %s ' % (f_ufahrej,f_ufahrej_screen)
else:
	cmd  = 'fselect %s %s ' % (f_ufahrej,f_ufahrej_screen)
	cmd += '"'
	for detid in exclude_detid_list:
		cmd += '(DET_ID != %s)' % (detid)
		if exclude_detid_list.index(detid) < len(exclude_detid_list)-1:
			cmd += ' && '
	cmd += '"'
print(cmd);os.system(cmd)

cmd  = 'gzip %s' % f_ufahrej_screen
print(cmd);os.system(cmd)
f_ufahrej_screen += '.gz'

# Make spectra and evaluate parameters for BG Mdeling
# 6.1 xselect: extract spectra, per GTI, from event lists
# for selected FPMs to obtain target spectra
# (from *cl.evt.gz or *cl50.evt.gz or selection-specialized files)
outdir_speclc = '%s/speclc' % outdir
cmd = 'rm -rf %s; mkdir -p %s' % (outdir_speclc,outdir_speclc)
print(cmd);os.system(cmd)

f_input_to_script = '%s/%s_gti_input.lst' % (outdir_speclc,prefix)
fin = open(f_input_to_script,'w')
f_gti_list2 = '%s/%s_gti.lst' % (outdir_speclc,prefix)
flog = open(f_gti_list2,'w')
dump = '# START STOP ID ObsID EXP(s) NFPM IBG(cps) HREJ(cps) NZ(cps) IBG52 HREJ52 NZ52'
flog.write(dump+'\n')
for line in open(f_gti_list):
	cols = line.split()
	tstart = float(cols[0])
	tstop  = float(cols[1])
	gtinum = int(cols[2])
	obsid = cols[3]
	texp   = float(cols[4])
	subdir = '%s/gti%05d' % (outdir_speclc,gtinum)
	cmd = 'rm -rf %s; mkdir -p %s' % (subdir,subdir)
	print(cmd);os.system(cmd)

	f_clevt_screen_gti_evt = '%s/%s_clmerge_screen_gti%d.evt' % (subdir,prefix,gtinum)
	f_clevt_screen_gti_pha = f_clevt_screen_gti_evt.replace('.evt','.pha')
	#cmd  = 'niextract-events '
	#cmd += 'filename="%s[TIME=%.1f:%.1f]" ' % (f_clevt_screen,tstart,tstop)
	#cmd += 'eventsout="%s" gti=GTI ' % f_clevt_screen_gti_evt
	#print(cmd);os.system(cmd)

	cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
set phaname PI
extract spectrum
save spectrum
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" % (f_clevt_screen,tstart,tstop,f_clevt_screen_gti_pha,subdir)
	print(cmd);os.system(cmd)

	f_ufatsel_screen_gti_ibg_lc = '%s/%s_ufatsel_screen_gti%d_ibg.lc' % (subdir,prefix,gtinum)
	cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 1500 1700
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" %  (f_ufatsel_screen,tstart,tstop,f_ufatsel_screen_gti_ibg_lc,subdir)
	print(cmd);os.system(cmd)
	gti_ibg_average = get_lc_rate(f_ufatsel_screen_gti_ibg_lc)

	f_ufahrej_screen_gti_hrej_lc = '%s/%s_ufahrej_screen_gti%d_hrej.lc' % (subdir,prefix,gtinum)
	cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 300 1800
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
mv xsel_timefile.asc xselect.log %s
""" %  (f_ufahrej_screen,tstart,tstop,f_ufahrej_screen_gti_hrej_lc,subdir)
	print(cmd);os.system(cmd)
	gti_hrej_average = get_lc_rate(f_ufahrej_screen_gti_hrej_lc)

	f_ufatsel_screen_gti_nz_lc = '%s/%s_ufatsel_screen_gti%d_nz.lc' % (subdir,prefix,gtinum)
	cmd = """
xselect <<EOF
xsel
read event %s .
yes
filter time scc 
%.4f,%.4f
x 
filter pha_cutoff 0 20
set binsize 1
show filter 
extract curve
save curve
%s
exit
no
EOF
rm -f xsel_timefile.asc xselect.log
mv xsel_timefile.asc xselect.log %s
""" %  (f_ufatsel_screen,tstart,tstop,f_ufatsel_screen_gti_nz_lc,subdir)
	print(cmd);os.system(cmd)
	gti_nz_average = get_lc_rate(f_ufatsel_screen_gti_nz_lc)

	dump = '%s %d %.4f %.4f %.4f' % (f_clevt_screen_gti_pha,num_of_fpm,
		gti_ibg_average,gti_hrej_average,gti_nz_average)
	fin.write(dump+'\n')

	dump  = '%.1f %.1f %d %s %.1f ' % (tstart,tstop,gtinum,obsid,texp)
	dump += '%d %.4f %.4f %.4f ' % (num_of_fpm,gti_ibg_average,gti_hrej_average,gti_nz_average)
	dump += '%.4f %.4f %.4f ' % (
		52.0 / float(num_of_fpm) * gti_ibg_average,
		52.0 / float(num_of_fpm) * gti_hrej_average,
		52.0 / float(num_of_fpm) * gti_nz_average
		)
	flog.write(dump+'\n')
fin.close()
flog.close()

expmap_IBGvsHREJ = ExposureMap_IBGvsHREJ()
expmap_NZ = ExposureMap_NZ()
f_gti_list3 = '%s/%s_gti.lst' % (outdir,prefix)
flog = open(f_gti_list3,'w')
dump = '# START STOP ID ObsID EXP(s) NFPM IBG(cps) HREJ(cps) NZ(cps) IBG52 HREJ52 NZ52 i_IBG j_HREJ k_NZ'
flog.write(dump+'\n')
for line in open(f_gti_list2):
	cols = line.split()
	if cols[0] == '#':
		continue 
	texp   = float(cols[4])
	NFPM = float(cols[5])	
	IBG52  = float(cols[9])
	HREJ52 = float(cols[10])
	NZ52   = float(cols[11])
	print(IBG52,HREJ52,NZ52,texp,NFPM)
	if not (IBG52 <= 10.0 and HREJ52 <= 5.0 and NZ52 <= 1600.0):
		print("skip...")
		continue
	i, j = expmap_IBGvsHREJ.add(IBG52,HREJ52,texp,NFPM)
	k = expmap_NZ.add(NZ52,texp,NFPM)
	flog.write(line.strip() + '%d %d %d\n' % (i, j, k))
expmap_IBGvsHREJ.show()
expmap_NZ.show()
flog.close()

norm_expmap_IBGvsHREJ = expmap_IBGvsHREJ.get_normalized_exposuremap() 
norm_expmap_NZ = expmap_NZ.get_normalized_exposuremap()
print(norm_expmap_IBGvsHREJ)
print(norm_expmap_NZ)

expr_list = []
for i in range(6):
	for j in range(5):
		bg_pha = '%s/bg_group_3C50_ngt_%d%d.pha' % (BGD_LIBRARY_PATH,i+1,j+1)
		print(i,j,norm_expmap_IBGvsHREJ[i][j],bg_pha)
		if norm_expmap_IBGvsHREJ[i][j] > 0.0:
			expr_list.append('%s * %.6f ' % (bg_pha,norm_expmap_IBGvsHREJ[i][j]))

for k in range(1,14):
	noise_pha = '%s/day_noise_nz%02d.pha' % (BGD_LIBRARY_PATH,k)
	print(k,norm_expmap_NZ[k],noise_pha)
	if norm_expmap_NZ[k] > 0.0:
		expr_list.append('%s * %.6f' % (noise_pha,norm_expmap_NZ[k]))
print(expr_list)		

expr = ''
for i in expr_list:
	expr += '%s' % i
	if expr_list.index(i) < len(expr_list)-1:
		expr += ' + '
print(expr)		

cmd  = 'mathpha '
cmd += '"%s" ' % expr
cmd += 'R outfil="test.pha" '
cmd += 'exposure=1.0 '
cmd += 'errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=yes'
print(cmd)

