#!/usr/bin/env python

__name__    = 'check_outlier_of_fpm_rate'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Apr. 16'

import os
import sys
import pyfits
import numpy as np
from optparse import OptionParser

NUM_OF_MPU = 7
NUM_OF_FPM = 8
NUM_OF_CUT_MIN = 1 
NUM_OF_CUT_MAX = 1 

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

help_message = """
(example)
check_outlier_of_fpm_rate.py -i 1012040144/xti/event_cl/ni1012040144_0mpu7_cl.evt.gz \
	-t 40 -d 0.5 -u 1.0 
"""
parser = OptionParser()
parser.add_option("-i","--inputevtfits",dest="inputevtfits",
	action="store",help="input event fits file",type="string")
parser.add_option("-o","--outputlc",dest="outputlc",#
	action="store",help="output lc file",type="string")
parser.add_option("-t","--timebin",dest="timebin",default=40.0,
	action="store",help="time bin size (sec)",type="float")
parser.add_option("-d","--eminkeV",dest="eminkeV",default=0.5,
	action="store",help="energy min (keV)",type="float")
parser.add_option("-u","--emaxkeV",dest="emaxkeV",default=1.0,
	action="store",help="energy max (keV)",type="float")
parser.add_option("-x","--thresholdsigma",dest="thresholdsigma",default=20,
	action="store",help="Threshold sigma",type="float")
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

pimin = int(nikeV2PI(options.eminkeV))
pimax = int(nikeV2PI(options.emaxkeV))


basename = os.path.splitext(os.path.basename(options.inputevtfits).replace('.gz',''))[0]

tmpdir = 'check_outlier_%s' % basename
cmd = 'rm -rf %s; mkdir -p %s' % (tmpdir, tmpdir)
print(cmd);os.system(cmd)

for mpu in range(0,NUM_OF_MPU):
	for fpm in range(0,NUM_OF_FPM):
		detid = "%02d" % (10*mpu+fpm)
		tmpoutfile = '%s/%s_detid%s.evt' % (tmpdir,basename,detid)
		cmd = 'fselect %s %s expr="(DET_ID==%s).and.(PI>=%d).and.(PI<%d)"' % (
			options.inputevtfits,tmpoutfile,detid,pimin,pimax)
		print(cmd);os.system(cmd)

		tmpoutlc = '%s/%s_detid%s_t%ds_%sto%skeV.flc' % (tmpdir,basename,detid,
			options.timebin,
			str(options.eminkeV).replace('.','p'),str(options.emaxkeV).replace('.','p'))
		cmd  = 'xselect <<EOF\n'
		cmd += 'xsel\n'
		cmd += 'read event %s .\n' % tmpoutfile
		cmd += 'yes\n'
		cmd += 'set binsize %d\n' % options.timebin
		cmd += 'extract curve\n'
		cmd += 'save curve\n'
		cmd += '%s\n' % tmpoutlc
		cmd += 'exit\n'
		cmd += 'no\n'
		cmd += 'EOF'
		print(cmd);os.system(cmd)
		cmd  = 'rm -f xselect.log'
		print(cmd);os.system(cmd)

		cmd  = 'ftcalc %s %s ' % (tmpoutlc,tmpoutlc)
		cmd += 'TIME "TIME+#TIMEZERO" ' 
		cmd += 'clobber=yes \n'
		cmd += 'ftcalc %s %s ' % (tmpoutlc,tmpoutlc)
		cmd += 'RATE%s RATE ' % detid 
		cmd += 'clobber=yes \n'		
		cmd += 'fdelcol %s+1 RATE N yes\n' % tmpoutlc		
		cmd += 'fdelcol %s+1 ERROR N yes\n' % tmpoutlc
		cmd += 'fdelcol %s+1 FRACEXP N yes\n' % tmpoutlc	
		print(cmd);os.system(cmd)

		if detid == "00":
			cmd  = 'cp %s %s ' % (tmpoutlc,options.outputlc)
			print(cmd);os.system(cmd)
		else:
			cmd  = 'fdelcol %s+1 TIME N yes\n' % tmpoutlc	
			cmd += 'ftpaste %s+1 %s+1 %s clobber=yes \n' % (options.outputlc,tmpoutlc,options.outputlc)
			print(cmd);os.system(cmd)

mean_list  = []
std_list   = []
maxv_list  = []
sigma_list = []

hdu = pyfits.open(options.outputlc)
for i in range(len(hdu['RATE'].data)):
	tline = hdu['RATE'].data[i]
	rate_detid_list = []
	for mpu in range(NUM_OF_MPU):
		for fpm in range(NUM_OF_FPM):
			rate_detid_list.append(float(tline['RATE%d%d' % (mpu,fpm)]))
	rate_detid_list.sort()
	cut_array = np.array(rate_detid_list[NUM_OF_CUT_MIN:-NUM_OF_CUT_MAX])
	mean = np.mean(cut_array)
	std  = np.std(cut_array)
	maxv = max(rate_detid_list)
	if mean == 0.0 and std == 0.0 and maxv == 0.0:
		sigma = 0.0
	elif std > 0.0:
		sigma = (maxv - mean)/std
	else:
		sigma = None

	print(i,mean,std,maxv,sigma)
	mean_list.append(mean)
	std_list.append(std)
	maxv_list.append(maxv)
	sigma_list.append(sigma)

cmd = 'rm -f tmp_sigma.fits'
print(cmd);os.system(cmd)	
cols = [
	pyfits.Column(name='MEAN',format='D',array=np.array(mean_list,'float')),
	pyfits.Column(name='STD',format='D',array=np.array(std_list,'float')),		
	pyfits.Column(name='MAXV',format='D',array=np.array(maxv_list,'float')),			
	pyfits.Column(name='SIGMA',format='D',array=np.array(sigma_list,'float'))
	]
hdu = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))
hdu.writeto('tmp_sigma.fits')

cmd = 'ftpaste %s+1 tmp_sigma.fits+1 %s clobber=yes \n' % (options.outputlc,options.outputlc)
print(cmd);os.system(cmd)
cmd = 'rm -rf tmp_sigma.fits %s' % tmpdir
print(cmd);os.system(cmd)

num_of_badtbin = np.sum(np.array(sigma_list,'float')<options.thresholdsigma)
cmd  = 'fparkey add=yes %d %s BADTBIN;' % (num_of_badtbin,options.outputlc)
cmd += 'fparkey add=yes %.3e %s BADSIGMA;' % (options.thresholdsigma,options.outputlc)
print(cmd);os.system(cmd)

cmd = 'rm -f gti_columns.txt gti_header.txt'
print(cmd);os.system(cmd)

f = open('gti_columns.txt','w')
dump = """START D s
STOP D s

"""
f.write(dump)
f.close()

f = open('gti_header.txt','w')
dump = """XTENSION= 'BINTABLE'           / binary table extension
BITPIX  =                    8 / 8-bit bytes
NAXIS   =                    2 / 2-dimensional binary table
NAXIS1  =                   16 / width of table in bytes
NAXIS2  =                    1 / number of rows in table
PCOUNT  =                    0 / size of special data area
GCOUNT  =                    1 / one data group (required keyword)
TFIELDS =                    2 / number of fields in each row
TTYPE1  = 'START   '           / lower GTI boundary
TFORM1  = 'D       '           / data format of field: 8-byte DOUBLE
TUNIT1  = 's       '           / physical unit of field
TTYPE2  = 'STOP    '           / upper GTI boundary
TFORM2  = 'D       '           / data format of field: 8-byte DOUBLE
TUNIT2  = 's       '           / physical unit of field
EXTNAME = 'STDGTI  '           / The name of this table
HDUCLASS= 'OGIP    '           / format conforms to OGIP standard
HDUCLAS1= 'GTI     '           / table contains Good Time Intervals
HDUCLAS2= 'STANDARD'           / standard Good Time Interval table
ONTIME  = 0.00000000000000E+00 / [s] sum of all Good Time Intervals
TSTART  = 0.00000000000000E+00 / [s] Lower bound of first GTI
TSTOP   = 0.00000000000000E+00 / [s] Uppler bound of last GTI
TIMEUNIT= 's       '           / All times in s unless specified otherwise
TIMESYS = 'TT      '           / XMM time will be TT (Terrestial Time)
TIMEREF = 'LOCAL   '           / Reference location of photon arrival times
TASSIGN = 'SATELLITE'          / Location of time assignment
TIMEZERO=                    0 / Clock correction (if not zero)
CLOCKAPP=                    T / Clock correction applied?
MJDREFI =                56658 / MJD reference day
MJDREFF = 7.775925925925930E-04 / MJD reference (fraction of day)
"""
f.write(dump)
f.close()

outbasename = os.path.splitext(os.path.basename(options.outputlc))[0]
if os.path.dirname(options.outputlc) in ["","/"]:
	bad_gti  = '%s_bad.gti' % (outbasename)
	good_gti = '%s_good.gti' % (outbasename)
else:
	bad_gti  = '%s/%s_bad.gti' % (os.path.dirname(options.outputlc),outbasename)
	good_gti = '%s/%s_good.gti' % (os.path.dirname(options.outputlc),outbasename)

cmd  = 'ftcopy "%s[1][SIGMA>%.3f]" tmp_bad.gti clobber=yes\n' % (options.outputlc,options.thresholdsigma)
cmd += 'ftcalc tmp_bad.gti tmp_bad.gti TSTART "TIME-0.5*(#TIMEDEL)" clobber=yes\n' 
cmd += 'ftcalc tmp_bad.gti tmp_bad.gti TSTOP "TIME+0.5*(#TIMEDEL)" clobber=yes\n' 
cmd += 'ftlist tmp_bad.gti+1 columns=TSTART,TSTOP rownum=no colheader=no opt=t > tmp_bad.txt'
print(cmd);os.system(cmd)

cmd  = 'ftcreate gti_columns.txt tmp_bad.txt %s headfile=gti_header.txt extname="GTI" clobber=yes\n' % bad_gti
cmd += 'rm -f tmp_bad.gti tmp_bad.txt\n'
print(cmd);os.system(cmd)

cmd  = 'ftcopy "%s[1][SIGMA<=%.3f]" tmp_good.gti clobber=yes\n' % (options.outputlc,options.thresholdsigma)
cmd += 'ftcalc tmp_good.gti tmp_good.gti TSTART "TIME-0.5*(#TIMEDEL)" clobber=yes\n' 
cmd += 'ftcalc tmp_good.gti tmp_good.gti TSTOP "TIME+0.5*(#TIMEDEL)" clobber=yes\n' 
cmd += 'ftlist tmp_good.gti+1 columns=TSTART,TSTOP rownum=no colheader=no opt=t > tmp_good.txt'
print(cmd);os.system(cmd)

cmd  = 'ftcreate gti_columns.txt tmp_good.txt %s headfile=gti_header.txt extname="GTI" clobber=yes\n' % good_gti
cmd += 'rm -f tmp_good.gti tmp_good.txt\n'
print(cmd);os.system(cmd)



cmd = 'rm -f gti_columns.txt gti_header.txt'
print(cmd);os.system(cmd)

