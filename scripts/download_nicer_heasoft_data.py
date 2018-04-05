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

parser = OptionParser()
parser.add_option("-c","--inputcsv",dest="inputcsv",
        action="store",help="Input NICER segment table (CSV format)",type="string")
parser.add_option("-k","--flag_kyoto",action="store_true",dest="flag_kyoto", default=True)
parser.add_option("-m","--flag_move",action="store_true",dest="flag_move", default=True)
(options, args) = parser.parse_args()

if len(args) != 1:
	sys.stderr.write('usage: %s.py obsid --inputcsv nicer_segment_table.csv\n' % __name__)	
	sys.stderr.write('   or: %s.py target_name --inputcsv nicer_segment_table.csv\n' % __name__)	
	quit()

if not os.path.exists(options.inputcsv):
	sys.stderr.write('inputcsv file does not exist.\n')
	quit()

df = pd.DataFrame.from_csv(options.inputcsv)

if args[0].isdigit():
	mode_obsid = True 
	target_obsid = str(args[0])
	source = df[df['Observation ID'] == int(target_obsid)]
	print("Target ObsID mode: {}".format(target_obsid))	
else:
	mode_obsid = False
	target_name = str(args[0])
	source = df[df['Target Name'] == target_name]	
	print("Target Name mode: {}".format(target_name))	
#print(source)


for n,[no, row] in enumerate(source.iterrows()):
	# Extract the relevant columns
	#obsid = "{:010d}".format(row['Observation ID'])
	obsid = str(row['Observation ID'])
	dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
	gexpo = row['Good Expo[s]']
	yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)
	srcname = row['Target Name']
	tarfile = '%s.tar' % obsid

	download_path = '%s/%s/%s ' % (FTPPATH_NICER_DATA, yyyy_mm, tarfile)
	dump = '%s %s  %s  %.1f (s) %s\n' % (srcname,obsid,dtime,gexpo, download_path)
	sys.stdout.write(dump)

	cmd  = 'curl -O %s' % download_path
	if options.flag_kyoto: 
		cmd += '-x ftp-proxy.kuins.net:8080 --proxy-user anonymous:'
	print(cmd);os.system(cmd)

	if options.flag_move:
		cmd = 'tar zxvf %s;' % tarfile 
		print(cmd);os.system(cmd)

		target_dir = '%s/%s' % (os.getenv('NICER_PUBLIC_DATA_PATH'),yyyy_mm)
		if not os.path.exists(target_dir):
			cmd = 'mkdir -p %s' % target_dir
			print(cmd);os.system(cmd)

		cmd = 'mv %s %s' % (obsid,target_dir)
		print(cmd);os.system(cmd)		





