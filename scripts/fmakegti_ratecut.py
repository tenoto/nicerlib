#!/usr/bin/env python
# 
# Original script was made by Slavko
"""
I think I've worked out a relatively straightforward way to filter a
NICER event list using a count rate cut from a binned light curve
produced from the same event list. This approach follows the same
procedure commonly employed for filtering out soft proton flares in
XMM-Newton data. The XMM SAS has command line tools to do this easily
but there are no equivalent FTOOLS, so I MacGyver'd a solution using a
combination of FTOOLS that has the same end result. This filtering
procedure can be useful for faint sources that are known to not be
strongly variable and should be more effective than making cuts based
just on cutoff rigidity or geographic location. For instance, during
the polar horn passages there may be occasional intervals of
relatively low background so, in principle, by using this method some
more usable data can be salvaged.

This is the basic procedure that I've found works without problems:
- First, generate a binned light curve, e.g. with Xselect. Based on it
decide what count rate threshold suits your science needs.
- Use ftcopy on the light curve to filter out time bins that have
count rates in excess of the threshold. For example, ftcopy
"lcurve_in.fits[1][RATE<2.0]" lcurve_out.fits will only keep the light
curve bins with count rates less than 2 counts/s
- With ftcalc compute the start and end times of the remaining time
bins, so they can be used as GTI intervals. (A side note here: time
zero in an Xselect light curve is actually the mid point of the first
bin and not the lower bound of the bin)
- Output the two computed table columns into an ascii file with ftlist
- Use ftcreate to generate a GTI FITS file using the ascii file of
start and end times, an ascii file listing the format of the data that
will go into the FITS file, and a bare-bones input header that is
required for the output FITS file to be recognized as a GTI file.
- The GTI file can then be applied to the NICER event list using
niextract-events by including the timefile='gti.fits[GTI]' keyword

I've put this sequence of FTOOL commands together in a bash script,
which is attached along with the template ascii files needed to
construct the GTI table file. As command line input it needs an event
file, a binned light curve, the time bin size of the light curve, the
desired count rate cut, and an output file name. Please test the
procedure out and feel free to hack the script for your own use. It
may be a good idea to include an implementation of it in Paul's
NICERsoft package.

# $1 - event file
# $2 - light curve file
# $3 - desired count rate cut (counts/s)
# $4 - time bin width
# $5 - output filtered event file

if [[ $# -ne 5 ]] ; then
    tput setaf 1 ; tput bold;   echo 'TOO FEW ARGUMENTS!'
    echo ' '
    tput sgr0; tput bold
    echo 'Usage: ./nicer_lc_filter.sh eventfile lcfile countrate timebin outfile'
    echo ' '
    tput sgr0
    echo 'Positional arguments:'
    echo '  eventfile          Input FITS event file name.'
    echo '  lcfile             Binned light curve input FITS file used for filtering' 
    echo '  countrate          Desired count rate cut (in counts/s). All light curve bins above this rate are removed.'
    echo '  timebin            Time bin width (in seconds) in input light curve.'
    echo '  outfile            Output filtered event list.'
    echo ' '
    tput sgr0
    exit 1 
fi

tput bold
echo "Making count rate cut of "$3" counts/s on events from "$1" using bins of width "$4" from light curve file "$2"."
tput sgr0
echo ' '

ftcopy "$2[1][RATE<$3]" lcurve_cut.fits clobber=yes

ftcalc lcurve_cut.fits lcurve_cut.fits TSTART "TIME-(0.5*$4)+#TIMEZERO" clobber=yes

ftcalc lcurve_cut.fits lcurve_cut.fits TEND "TIME+(0.5*$4)+#TIMEZERO" clobber=yes

ftlist lcurve_cut.fits[1] columns=TSTART,TEND rownum=no colheader=no opt=t > gti_data.txt

ftcreate gti_columns.txt gti_data.txt gti.fits headfile=gti_header.txt extname="GTI" clobber=yes

niextract-events $1 $5 timefile='gti.fits[GTI]'

echo "Filtered event lists written in "$5"."
"""

__name__    = 'fmakegti_ratecut'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Feb 20'

import os
import pyfits
from optparse import OptionParser
from datetime import datetime 

# $3 - desired count rate cut (counts/s)
# $4 - time bin width

parser = OptionParser()
parser.add_option("-i","--inputlcfits",dest="inputlcfits",
	action="store",help="input light curve fits file",type="string")
parser.add_option("-o","--outputgti",dest="outputgti",
	action="store",help="output gti file name",type="string")
parser.add_option("-d","--ratemin",dest="ratemin",
	action="store",help="Threshold rate minimum (cps)",type="float")
parser.add_option("-u","--ratemax",dest="ratemax",
	action="store",help="Threshold rate maximum (cps)",type="float")
#parser.add_option("-b","--timebin",dest="timebin",
#	action="store",help="Time bin width",type="float")
(options, args) = parser.parse_args()

if options.inputlcfits == None:
	print "input light curve fits file is needed. %> -i inputlcfits"
	quit()
if options.outputgti == None:	
	print "output gti file name is needed. %> -o outputgti"
	quit()
if options.ratemin == None:
	print "ratemin (cps) is needed. %> -d ratemin"
	quit()
if options.ratemax == None:
	print "ratemax (cps) is needed. %> -u ratemax"
	quit()	
#if options.timebin == None:
#	print "timebin (s) is needed. %> -b timebin"
#	quit()
print "inputlcfits : %s " % options.inputlcfits
print "outputgti: %s " % options.outputgti
print "ratemin: %.3e (cps)" % options.ratemin
print "ratemax: %.3e (cps)" % options.ratemax
#print "timebin: %.2f (s)" % options.timebin

if not os.path.exists(options.inputlcfits):
	print "input light curve fits file does not exists: %s" % options.inputlcfits
	quit()
if os.path.exists(options.outputgti):
	print "output gti file has already existed: %s " % options.outputgti
	quit()

hdu = pyfits.open(options.inputlcfits)
timebin = hdu[0].header['TIMEDEL']

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

cmd  = 'rm -f tmp_lcurve_cut.fits tmp_gti_data.txt tmp_gti_data_shrink.txt\n'
cmd += 'ftcopy "%s[1][(RATE>=%.3e)&&(RATE<=%.3e)]" tmp_lcurve_cut.fits clobber=yes\n' % (options.inputlcfits,options.ratemin,options.ratemax)
print(cmd);os.system(cmd)
hdu = pyfits.open('tmp_lcurve_cut.fits')

f = open('tmp_gti_data_shrink.txt','w')
if len(hdu['RATE'].data) <= 0:
	print("--- no GTI...")
	dump = '0.0 0.0\n'
	f.write(dump)
else:
	print("--- GTI exists.")	
	cmd  = 'ftcalc tmp_lcurve_cut.fits tmp_lcurve_cut.fits TSTART "TIME-(0.5*%.3f)+#TIMEZERO" clobber=yes\n' % (timebin)
	cmd += 'ftcalc tmp_lcurve_cut.fits tmp_lcurve_cut.fits TEND "TIME+(0.5*%.3f)+#TIMEZERO" clobber=yes\n' % (timebin)
	cmd += 'ftlist tmp_lcurve_cut.fits[1] columns=TSTART,TEND rownum=no colheader=no opt=t > tmp_gti_data.txt\n'
	print(cmd);os.system(cmd)
	flag_first = True
	for line in open('tmp_gti_data.txt'):
		cols = line.split()
		tmp_TSTART = cols[0]
		tmp_TSTOP  = cols[1]
		if flag_first:
			TSTART = tmp_TSTART
			prev_TSTOP  = tmp_TSTOP
			flag_first = False
			continue

		if tmp_TSTART != prev_TSTOP:
			dump = '%s %s\n' % (TSTART,prev_TSTOP)
			f.write(dump)
			TSTART = tmp_TSTART 
		prev_TSTOP = tmp_TSTOP
f.close()

cmd  = 'ftcreate gti_columns.txt tmp_gti_data_shrink.txt %s headfile=gti_header.txt extname="GTI" clobber=yes\n' % options.outputgti
cmd += 'rm -f tmp_lcurve_cut.fits tmp_gti_data.txt tmp_gti_data_shrink.txt gti_columns.txt gti_header.txt\n'
cmd += 'fparkey %.5e %s RATEMIN add=yes comm="Rate minimum (cps)"\n' % (options.ratemin,options.outputgti)
cmd += 'fparkey %.5e %s RATEMAX add=yes comm="Rate maximum (cps)"\n' % (options.ratemax,options.outputgti)
print(cmd);os.system(cmd)

