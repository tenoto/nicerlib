#!/usr/bin/env python

#wget --passive-ftp -q ftp://heasarc.gsfc.nasa.gov/nicer/data/obs/2017_08/1020270115.tar 

__name__    = 'download_nicer_heasoft_data'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 5'

FTPPATH_NICER_DATA = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs'
DELAY_DAYS_FOR_PUBLIC = 14

import os 
import sys 
from optparse import OptionParser
import pandas as pd
import datetime 
import subprocess

parser = OptionParser()
parser.add_option("-c","--inputcsv",dest="inputcsv",
        action="store",help="Input NICER segment table (CSV format)",type="string")
parser.add_option("-s","--skip_yyyymm_list",dest="skip_yyyymm_list",default=None,
        action="store",help="skip list of [yyyymm], e.g., \"[2018_02,2018_03]\"",type="string")
parser.add_option("-m","--flag_move",action="store_true",dest="flag_move",default=True)
(options, args) = parser.parse_args()

if len(args) != 1:
	sys.stderr.write('usage: %s.py obsid --inputcsv nicer_segment_table.csv\n' % __name__)	
	sys.stderr.write('   or: %s.py target_name --inputcsv nicer_segment_table.csv\n' % __name__)	
	quit()

if not os.path.exists(options.inputcsv):
	sys.stderr.write('inputcsv file does not exist.\n')
	quit()

df = pd.DataFrame.from_csv(options.inputcsv)
df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)

if args[0].isdigit():
	mode_obsid = True 
	target_obsid = str(args[0])
	source = df[df['Observation ID'] == str(target_obsid)]
	print("Target ObsID mode: {}".format(target_obsid))	
else:
	mode_obsid = False
	target_name = str(args[0])
	source = df[df['Target Name'] == target_name]	
	print("Target Name mode: {}".format(target_name))	
#print(source)

if options.skip_yyyymm_list is not None:
	skip_yyyymm_list = list(options.skip_yyyymm_list.replace('[','').replace(']','').split(','))

for n,[no, row] in enumerate(source.iterrows()):
	# Extract the relevant columns
	#obsid = "{:010d}".format(row['Observation ID'])
	obsid = str(row['Observation ID'])
	try:
		dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
	except:
		print("... skip (ObsID=%s) because Start TimeUTC is blank." % obsid)
		continue 
	gexpo = row['Good Expo[s]']
	yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)
	srcname = row['Target Name']
	tarfile = '%s.tar' % obsid

	if options.skip_yyyymm_list is not None:
		if yyyy_mm in skip_yyyymm_list:
			print("... skip yyyy_mm.")
			continue

	obsid_type = obsid[0]
	if obsid_type == '0':
		print("... skip ObsID starting 0 (ObsID=%s)." % obsid)
		continue 

	t_now = datetime.datetime.now()
	if (t_now - dtime).days < DELAY_DAYS_FOR_PUBLIC:
		print("... skip (ObsID=%s) since the observation is within %d days" % (obsid,DELAY_DAYS_FOR_PUBLIC))
		continue 

	target_dir = '%s/%s' % (os.getenv('NICER_PUBLIC_DATA_PATH'),yyyy_mm)

	download_path = '%s/%s/%s' % (FTPPATH_NICER_DATA, yyyy_mm, obsid)
	dump = '%s %s  %s  %.1f (s) %s\n' % (srcname,obsid,dtime,gexpo, download_path)
	sys.stdout.write(dump)

	for subdir in ['auxil','log','xti']:
		cmd  = 'wget -q -nH --no-check-certificate --cut-dirs=5 '
		cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
		cmd += '%s/%s' % (download_path,subdir)
		print(cmd);os.system(cmd)

	if options.flag_move:
		if not os.path.exists(target_dir):
			cmd = 'mkdir -p %s' % target_dir
			print(cmd);os.system(cmd)

		cmd = 'mv %s %s' % (obsid,target_dir)
		print(cmd);os.system(cmd)		




