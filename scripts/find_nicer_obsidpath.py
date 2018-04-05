#!/usr/bin/env python

__name__    = 'find_nicer_obsid'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 5'

import os 
import sys 
from optparse import OptionParser
import pandas as pd
import datetime 

usage = "usage: %prog [options] target_name"
parser = OptionParser(usage=usage)
parser.add_option("-c","--inputcsv",dest="inputcsv",
        action="store",help="Input NICER segment table (CSV format)",type="string")
parser.add_option("-v","--verbose",action="store_true",dest="verbose", default=False)
(options, args) = parser.parse_args()

if not os.path.exists(options.inputcsv):
	sys.stderr.write('inputcsv file does not exist.\n')
	quit()
df = pd.DataFrame.from_csv(options.inputcsv)
df['Observation ID'] = df['Observation ID'].astype(str).str.zfill(10)

if len(args) != 1:
	sys.stderr.write('usage: %s.py target_name --inputcsv nicer_segment_table.csv\n' % __name__)	
	quit()

target_name = args[0]

selected = df[df['Target Name']==target_name]
#target_id = selected.head(1)['Target ID']
#print("==============================================================")
#print("Target Name %s (Target ID %d)" % (target_name, target_id))
#print("--------------------------------------------------------------")	
#if options.verbose:
#	print(selected.loc[:,['Start TimeUTC','Observation ID','Stop TimeUTC','On-Targ Expo[s]','Good Expo[s]','Process State','Process Date']])
#else:
#	print(selected.loc[:,['Start TimeUTC','Observation ID','On-Targ Expo[s]','Good Expo[s]']])
#print("--------------------------------------------------------------")	
#print("Target Name %s (Target ID %d)" % (target_name, target_id))
#print("==============================================================")

for n,[no, row] in enumerate(selected.iterrows()):
	obsid  = str(row['Observation ID'])
	dtime = datetime.datetime.strptime(row['Start TimeUTC'], "%Y-%m-%dT%H:%M:%S")
	yyyy_mm = "{}_{:02d}".format(dtime.year, dtime.month)	
	obsid_path = '%s/%s/%s' % (os.getenv('NICER_PUBLIC_DATA_PATH'),yyyy_mm,obsid)
	if os.path.exists(obsid_path):
		print(obsid_path)
