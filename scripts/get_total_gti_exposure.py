#!/usr/bin/env python

import sys 
import pyfits

if len(sys.argv) != 2:
	sys.stderr.write("%> %s" % sys.argv[0])
	quit()
infits = sys.argv[1]

total_exp = 0.0
hdu = pyfits.open(infits)
for line in hdu['GTI'].data:
	exp = line['STOP'] - line['START']
	#print(line['START'],line['STOP'],exp)
	total_exp += exp
print(total_exp)
