#!/usr/bin/env python

__name__    = 'fxselect_filter_energy'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 March 30'

import os
import sys
#import pyfits
import astropy.io.fits as pyfits
from optparse import OptionParser
from datetime import datetime 

# Conversion from PI to keV (PI is in units of 10 eV)
PI_TO_KEV = 10.0/1000.0
KEV_TO_PI = 1000.0/10.0

parser = OptionParser()
parser.add_option("-i","--inputevtfits",dest="inputevtfits",
	action="store",help="input event fits file",type="string")
parser.add_option("-o","--outputevtfits",dest="outputevtfits",
	action="store",help="output lc file",type="string")
parser.add_option("-d","--emin",dest="emin",
	action="store",help="Energy min (keV)",type="float")
parser.add_option("-u","--emax",dest="emax",
	action="store",help="Energy max (keV)",type="float")
(options, args) = parser.parse_args()

if options.inputevtfits == None:
	print "input event fits file is needed. %> -i inputevtfits"
	quit()
if options.outputevtfits == None:	
	print "output lc file is needed. %> -o outputevtfits"
	quit()
if options.emin == None:	
	print "emin (keV) is needed. %> -u emin"
	quit()
if options.emax == None:	
	print "emax (keV) is needed. %> -u emax"
	quit()

print "input event fits file  : %s " % options.inputevtfits
print "output event fits file : %s " % options.outputevtfits	
print "emin (keV) : %.2f " % options.emin	
print "emax (keV) : %.2f " % options.emax

if not os.path.exists(options.inputevtfits):
	print "input event fits file does not exists: %s" % options.inputevtfits
	quit()
if os.path.exists(options.outputevtfits):
	print "output event fits file has already existed: %s " % options.outputevtfits
	quit()

hdu = pyfits.open(options.inputevtfits)
if len(hdu['EVENTS'].data) == 0:
	sys.stdout.write('Skip: No events in fits file %s' % options.inputevtfits)
	quit()

pi_min = int(KEV_TO_PI * options.emin)
pi_max = int(KEV_TO_PI * options.emax)
cmd  = 'fselect %s %s ' % (options.inputevtfits,options.outputevtfits)
cmd += '"(PI >= %d).and.(PI < %d)" ' % (pi_min,pi_max)
print(cmd);os.system(cmd)
#return options.outputevtfits
