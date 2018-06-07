#!/usr/bin/env python

__name__    = 'generate_mitbgd_3C50'
__author__  = 'Teru Enoto (for Rons model)'
__version__ = '1.00'
__date__    = '2018 June 6'

"""
see the document 'NICER XTI Background Model 3C50'
Description and User Instructions
(BGMod_3C50.tar)
v.3 2018 0531
"""

BGD_LIBRARY_PATH = "/Users/enoto/Dropbox/enoto/research/nicer/data/mitbgd/BGMod_3C50"
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
parser.add_argument('--checkoutlier', action='store_true', 
	help='flag to check outlier.') 
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

obsid_list = []
ufaevt_list = []
clevt_list = []
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
	hdu = pyfits.open(clevt)
	print(clevt)
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
				gti_list.append([nid,tstart,tstop,tgti])
				message = '  %d %.1f %.1f %.1f ' % (nid,tstart,tstop,tgti)
				f.write('%.1f %.1f\n' % (tstart,tstop))
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
for gti in gti_list:
	gtinum = gti[0]
	tstart = gti[1]
	tstop  = gti[2]
	texp   = gti[3]
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

fin.close()


"""
		orgdir = os.getcwd()
		os.chdir(self.bgddir)
		modeldir = '%s/BGMod_3C50' % os.getenv('NICER_BGDMODEL_PATH')
		syslink_files = ['bg_group_3C50.table','day_noise_3C50.table','bg_group_3C50_ngt_*.pha','day_noise_nz*.pha','make_BG3C50_spectrum.go']
		for file in syslink_files:
			cmd = 'ln -s %s/%s .' % (modeldir,file)
			print(cmd);os.system(cmd)

		cmd  = 'rm -f setup_BG3C50.go;'
		cmd += 'awk \'{print "./make_BG3C50_spectrum.go", $1, $2, $3, $4, $5 " ;'
		cmd += 'cat calc_bg.go >> calc_BG3C50.xcm"}\' %s > setup_BG3C50.go' % os.path.basename(self.inputfile_bgd_param)
		print(cmd);os.system(cmd)

		cmd  = 'source setup_BG3C50.go'
		print(cmd);os.system(cmd)

		cmd  = 'source calc_BG3C50.xcm'
		print(cmd);os.system(cmd)

		for srcpha in glob.glob('ni%s_gti*.pha' % self.obsid):
			hdu = pyfits.open(srcpha)
			exposure = float(hdu['SPECTRUM'].header['EXPOSURE'])
			bgdpha = 'bg_%s' % srcpha
			#for extnum in range(3):
			for extnum in [1]:
				cmd = 'fparkey %.7f %s+%d EXPOSURE' % (exposure,bgdpha,extnum)
				print(cmd);os.system(cmd)

		for file in syslink_files:
			cmd = 'rm -f %s' % (file)
			print(cmd);os.system(cmd)			
		tentative_files = ['calc_bg.go','find_bg3C50_lib_spec.go','find_group.go','find_nzgroup.go','group.dat','inp.norms','nzgroup.dat','temp*','test*']
		for file in tentative_files:
			cmd = 'rm -f %s' % (file)
			print(cmd);os.system(cmd)			
		os.chdir(orgdir)
"""
