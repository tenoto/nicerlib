#!/usr/bin/env python

#wget --passive-ftp -q ftp://heasarc.gsfc.nasa.gov/nicer/data/obs/2017_08/1020270115.tar 

__name__    = 'download_nicer_heasoft_data'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 5'

FTPPATH_NICER_DATA = 'ftp://heasarc.gsfc.nasa.gov/nicer/data/obs'

import os 
import sys 
from optparse import OptionParser
import pandas as pd
import datetime 
import subprocess

parser = OptionParser()
parser.add_option("-c","--inputcsv",dest="inputcsv",
        action="store",help="Input NICER segment table (CSV format)",type="string")
parser.add_option("-k","--flag_kyoto",action="store_true",dest="flag_kyoto", default=True)
parser.add_option("-m","--flag_move",action="store_true",dest="flag_move", default=True)
parser.add_option("-s","--skip_yyyymm_list",dest="skip_yyyymm_list",default=None,
        action="store",help="skip list of [yyyymm], e.g., \"[2018_02,2018_03]\"",type="string")
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

skip_yyyymm_list = list(options.skip_yyyymm_list.replace('[','').replace(']','').split(','))

for n,[no, row] in enumerate(source.iterrows()):
	# Extract the relevant columns
	#obsid = "{:010d}".format(row['Observation ID'])
	obsid = str(row['Observation ID'])
	dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
	gexpo = row['Good Expo[s]']
	yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)
	srcname = row['Target Name']
	tarfile = '%s.tar' % obsid

	if yyyy_mm in skip_yyyymm_list:
		print("... skip yyyy_mm.")
		continue

	download_path = '%s/%s/%s ' % (FTPPATH_NICER_DATA, yyyy_mm, tarfile)
	dump = '%s %s  %s  %.1f (s) %s\n' % (srcname,obsid,dtime,gexpo, download_path)
	sys.stdout.write(dump)

	target_dir = '%s/%s' % (os.getenv('NICER_PUBLIC_DATA_PATH'),yyyy_mm)
	cmd  = 'curl -f -O %s' % download_path
	if options.flag_kyoto: 
		cmd += '-x ftp-proxy.kuins.net:8080 --proxy-user anonymous:'
	if os.path.exists('%s/%s' % (target_dir,obsid)):
		print("...skip, direcotry already exists...")
		continue 
	print(cmd);os.system(cmd)
	#resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	#print(resp.stderr.readlines())

	if options.flag_move:
		try:
			cmd = 'tar zxvf %s;' % tarfile 
			print(cmd);os.system(cmd)

			cmd = 'rm -f %s' % tarfile
			print(cmd);os.system(cmd)
			
			if not os.path.exists(target_dir):
				cmd = 'mkdir -p %s' % target_dir
				print(cmd);os.system(cmd)

			cmd = 'mv %s %s' % (obsid,target_dir)
			print(cmd);os.system(cmd)		
		except:
			print("error of opening the tarfile (download error).")





