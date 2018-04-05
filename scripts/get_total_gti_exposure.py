#!/usr/bin/env python

import sys 
import pyfits

if len(sys.argv) != 2:
	sys.stderr.write("%> %s" % sys.argv[0])
	quit()
infits = sys.argv[1]

total_exp = 0.0
hdu = pyfits.open(infits)
extname_list = []
for i in range(len(hdu)):
	try:
		extname_list.append(hdu[i].header['EXTNAME'])
	except:
		print("no EXTNAME in the extension %d" % i)

if 'GTI' in extname_list:
	gtiname = 'GTI'
elif 'STDGTI' in extname_list:
	gtiname = 'STDGTI'

for line in hdu[gtiname].data:
	exp = line['STOP'] - line['START']
	#print(line['START'],line['STOP'],exp)
	total_exp += exp
print(total_exp)
