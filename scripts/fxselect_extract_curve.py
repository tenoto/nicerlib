#!/usr/bin/env python

__name__    = 'fxselect_extract_curve'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Feb 21'

import os
import sys
import astropy.io.fits as pyfits
from optparse import OptionParser
from datetime import datetime 


def nikeV2PI(keV):
	"""
	# PI = ENERGY / 10 eV = ENERGY / 0.010 keV = 100 * ENERGY / keV 
	# ENERGY = 10 eV * PI = 0.010 keV * PI
	"""	
	return keV * 100.0 

def niPI2keV(PI):
	"""
	# PI = ENERGY / 10 eV = ENERGY / 0.010 keV = 100 * ENERGY / keV 
	# ENERGY = 10 eV * PI = 0.010 keV * PI
	"""		
	return PI / 100.0 

parser = OptionParser()
parser.add_option("-i","--inputevtfits",dest="inputevtfits",
	action="store",help="input event fits file",type="string")
parser.add_option("-o","--outputlc",dest="outputlc",
	action="store",help="output lc file",type="string")
parser.add_option("-t","--timebin",dest="timebin",default=10.0,
	action="store",help="time bin size (sec)",type="float")
parser.add_option("-d","--eminkeV",dest="eminkeV",default=0,
	action="store",help="energy min (keV)",type="float")
parser.add_option("-u","--emaxkeV",dest="emaxkeV",default=15.0,
	action="store",help="energy max (keV)",type="float")
(options, args) = parser.parse_args()

if options.inputevtfits == None:
	print "input event fits file is needed. %> -i inputevtfits"
	quit()
if options.outputlc == None:	
	print "output lc file is needed. %> -o outputlc"
	quit()

print "input event fits file  : %s " % options.inputevtfits
print "output lc file: %s " % options.outputlc	

if not os.path.exists(options.inputevtfits):
	print "input event fits file does not exists: %s" % options.inputevtfits
	quit()
if os.path.exists(options.outputlc):
	print "output lc file has already existed: %s " % options.outputlc
	quit()

hdu = pyfits.open(options.inputevtfits)
if len(hdu['EVENTS'].data) == 0:
	sys.stdout.write('Skip: No events in fits file %s' % options.inputevtfits)
	quit()

pimin = int(nikeV2PI(options.eminkeV))
pimax = int(nikeV2PI(options.emaxkeV))

cmd = """
xselect <<EOF
xsel
read event %s .
yes
set binsize %f
filter pha_cutoff %d %d 
show filter 
extract curve
save curve
%s
exit
no
EOF
""" % (options.inputevtfits,options.timebin,pimin,pimax,options.outputlc)
print(cmd);os.system(cmd)	

cmd = 'rm -f xselect.log'
print(cmd);os.system(cmd)

