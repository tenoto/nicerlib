#!/usr/bin/env python

__name__    = 'nitimeconv'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 6'

from optparse import OptionParser
from jdcal import MJD_0, MJD_JD2000, jd2gcal, gcal2jd
import datetime

MJDREFI   = 56658.0
MJDREFF   = 0.000777592592592593
TIMEZERO  = 0.0
LEAP_INIT = 2.0

usage = """ %prog -t inputtype -i value

NOTE:
* See the definition: https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
* NICER Time is reported in the TIME column.
	- TIME is elapsed TT seconds since the epoch 2014-01-01T00:00:00 UTC
* TSTART and TSTOP report the start and stop of good time 
* EXPOSURE is total live time 
* Conversion of NICER timestamps to absolute time in MJD can be done with the following:
	- MJD(TT) = MJDREFI+MJDREFF+(TIMEZERO+TIME)/ 86400
	- MJD(UTC) = MJD(TT) + LEAP_INIT
* BUT: See the Caveats pages for important limita0ons! Currently LEAP_INIT is not properly filled so you will need to use LEAP_INIT=2 manually.	
"""
parser = OptionParser(usage=usage)
parser.add_option("-t","--type",dest="type",default="mission",
       action="store",help="Any of mission/date",type="string")
parser.add_option("-i","--input",dest="input",
        action="store",help="Input time (sec)",type="string")
(options, args) = parser.parse_args()

mission = None
mjd_tt  = None
mjd_utc = None
jd_tt   = None
jd_utc  = None
date_tt = None

if options.type == "mission":
	mission = float(options.input)

	mjd_tt  = MJDREFI+MJDREFF+(TIMEZERO+mission)/86400.0
	mjd_utc = mjd_tt + LEAP_INIT
	jd_tt  = mjd_tt + MJD_0
	jd_utc = mjd_utc + MJD_0 

	yyyy, mm, dd, dayflac = jd2gcal(MJD_0,mjd_tt)
	tsec = dayflac* 24.0 * 60.0 * 60.0
	date_tt = '%s-%02d-%02d %s' % (yyyy, mm, dd, str(datetime.timedelta(seconds=tsec)))

elif options.type == "date_tt":
	date_tt = str(options.input)	

dump  = "MET    : %.8f (NICER TIME, Misson Elapsed Time)\n" % mission
dump += "MJD_TT : %.8f\n" % mjd_tt
dump += "MJD_UTC: %.8f\n" % mjd_utc
dump += "JD_TT  : %.8f\n" % jd_tt
dump += "JD_UTC : %.8f\n" % jd_utc
dump += "DATE_TT: %s\n" % date_tt
#dump += "DATE_TT: %s\n" % gcal_tt
print(dump)





