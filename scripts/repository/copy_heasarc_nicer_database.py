#!/usr/bin/env python

__name__    = 'copy_heasarc_nicer_database'
__author__  = 'Teru Enoto'
__version__ = '1.00'
__date__    = '2018 April 5'

import os 
import sys 
import glob 
import datetime 
import pandas as pd
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-u","--username",dest="username",
        action="store",help="NICER team username to access the webpage.",type="string")
parser.add_option("-p","--password",dest="password",
        action="store",help="NICER team password to access the webpage.",type="string")
parser.add_option("-n", "--newtable",action="store_true",
	dest="flag_download_newtable", default=False,
	help="flag to download a new segment table.")
(options, args) = parser.parse_args()

REMOTE_FTPPATH_NICER_DATA = 'https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs'
DIRPATH_LOCAL_NICER_DATA  = 'FTP/nicer/data/obs'

dir_input = 'script/input'
fname_original_download_source_list = '%s/scripts/repository/download_source_list.txt' % os.getenv('NICER_SOFT_PATH')
fname_logfile = '%s/download.log' % dir_input

#======================================
# Download segment 
#======================================
if not os.path.exists(dir_input):
	cmd = 'mkdir -p %s' % (dir_input)
	print(cmd);os.system(cmd)

if options.flag_download_newtable:
	date = datetime.datetime.now()
	datestr = date.strftime('%Y%m%d_%H%M')
	fname_segment_table = '%s/nicer_segment_table_v%s.csv' % (dir_input,datestr)

	print("...downloading the latest segment table...")
	cmd = 'rm -f %s/nicer_segment_table_v*.csv\n' % dir_input
	cmd += 'download_nicer_segment_table_to_csv.py '
	cmd += '--username %s ' % options.username
	cmd += '--password %s ' % options.password 
	cmd += '--outcsv %s ' % fname_segment_table
	print(cmd);os.system(cmd)
else:
	csvlst = glob.glob('%s/nicer_segment_table_v*.csv' % dir_input)
	if len(csvlst) == 0:
		sys.stderr.write('no nicer_segment_table_v*.csv\n')
		quit()
	else:
		fname_segment_table = csvlst[0]
print('fname_segment_table:%s' % fname_segment_table)

#======================================
# Prepare Download Source List
#======================================
cmd  = 'rm -f download_source_list.txt;\n'
cmd += 'cp %s %s' % (fname_original_download_source_list,dir_input)
print(cmd);os.system(cmd)

DOWNLOAD_SOURCE_LIST = []
for line in open(fname_original_download_source_list):
	cols = line.split()
	if cols[0] == '#':
		continue
	DOWNLOAD_SOURCE_LIST.append(cols[0])
print(DOWNLOAD_SOURCE_LIST)	

#======================================
# Prepare directories
#======================================
if not os.path.exists(DIRPATH_LOCAL_NICER_DATA):
	cmd = 'mkdir -p %s' % DIRPATH_LOCAL_NICER_DATA
	print(cmd);os.system(cmd)
if not os.path.exists('data'):
	cmd = 'ln -s FTP data'
	print(cmd);os.system(cmd)

#======================================
# Prepare Download Source List
#======================================
df = pd.DataFrame.from_csv(fname_segment_table)
df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)

f_logfile = open(fname_logfile,'a')
for index, row in df.iterrows():
	if not row['Target Name'] in DOWNLOAD_SOURCE_LIST:
		continue 
	if not row['Process State'] in ['DELV','DELVV']:
		continue 
	if row['Observation ID'][0] == '0':
		continue

	try:
		dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
	except:
		print("... skip (ObsID=%s) because Start TimeUTC is blank." % obsid)
		continue 
	yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)
	obsid = str(row['Observation ID'])
	tarfile = '%s.tar' % obsid
	target_path = '%s/%s/%s' % (REMOTE_FTPPATH_NICER_DATA, yyyy_mm, obsid)
	dir_move_to = '%s/%s' % (DIRPATH_LOCAL_NICER_DATA, yyyy_mm)

	dir_final_path = '%s/%s' % (dir_move_to,obsid)
	if os.path.exists(dir_final_path):
		sys.stderr.write('target obsid %s is already existed. skipped.\n' % obsid)
		continue

	if not os.path.exists(dir_move_to):
		cmd = 'mkdir -p %s' % dir_move_to
		print(cmd);os.system(cmd)

	now = datetime.datetime.now().isoformat()
	message = '%s: %s %s %s\n' % (now,row['Target Name'],row['Observation ID'],row['Start TimeUTC'])
	sys.stdout.write(message)
	f_logfile.write(message)
	try:
		#cmd = ''
		#for subdir in ['auxil','log','xti']:
		#	cmd += 'wget -q -nH --no-check-certificate --cut-dirs=5 '
		#	cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
		#	cmd += '%s/%s;' % (target_path,subdir)
		subdir = 'xti'
		cmd  = 'wget -q -nH --no-check-certificate --cut-dirs=5 '
		cmd += '-r -l0 -c -N -np -R \'index*\' -erobots=off --retr-symlinks '
		cmd += '%s/%s;' % (target_path,subdir)		
		cmd += 'mv %s %s;\n' % (obsid, dir_move_to)
		print(cmd);os.system(cmd)
		f_logfile.write(cmd)
		message = '...downloaded %s\n' % obsid
	except:
		message = '...error: download faield %s\n' % obsid
		sys.stderr.write(message)
	sys.stdout.write(message)
	f_logfile.write(message)
	#exit()
f_logfile.close()



