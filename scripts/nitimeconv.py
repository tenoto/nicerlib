#!/usr/bin/env python

__name__    = 'nitimeconv'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 6'

from optparse import OptionParser
from astropy.time import Time 

NICER_MJDREFI   = 56658.0
NICER_MJDREFF   = 0.000777592592592593
NICER_TIMEZERO  = 0.0
LEAP_INIT = 2.0

FERMI_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
NUSTART_MET_ORIGIN = Time('2010-01-01T00:00:00.000',format='isot',scale='utc')
RXTE_MET_ORIGIN    = Time('1994-01-01T00:00:00.000',format='isot',scale='utc')
SUZAKU_MET_ORIGIN  = Time('2000-01-01T00:00:00.000',format='isot',scale='utc')
SWIFT_MET_ORIGIN   = Time('2001-01-01T00:00:00.000',format='isot',scale='utc')
XMM_MET_ORIGIN     = Time('1998-01-01T00:00:00.000',format='isot',scale='tt')
CHANDRA_MET_ORIGIN = XMM_MET_ORIGIN
"""
Fermi seconds since 2001.0 UTC (decimal)		410227203.000
Fermi mission week (integer)		291
LIGO/GPS seconds since 1980-01-06 UTC (decimal)		1072569616.000
NuSTAR seconds since 2010.0 UTC (decimal)		126230401.000
RXTE seconds since 1994.0 UTC (decimal)		631152007.000
RXTE seconds since 1994.0 UTC (hexadecimal)		0x259e9d87
RXTE mission day number (dddd:hh:mm:ss)		7305:00:00:07.000
RXTE decimal mission day (dddd.ddd...)		7305.00008102
Suzaku seconds since 2000.0 UTC (decimal)		441849603.000
Swift seconds since 2001.0 UTC (decimal)		410227211.463
XMM/Chandra seconds since 1998.0 TT (decimal)		504921667.184
"""

usage = """ 
NAME 
	nitimeconv - Convert mission time in different time format and scale 

USAGE 
	%prog intime -f format -s scale 

DESCRIPTION
	'%prog' takes an input time with time format and scale (see below), 
	and converts the tine in another time systems. This gives similar outputs 
	of HEASARC "xTime - A Date/Time Conversion Utility" 
		https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
	This internal conversion is performed by pytho astropy library.
		http://docs.astropy.org/en/stable/api/astropy.time.Time.html#astropy.time.Time
	NICER Mission Elapse Time (MET) in the Time Column is defined as elapsed TT
	seconds since the epoch 2014-01-01T00:00:00 UTC (corresponding with 56658.00077759
	MJD_TT). Conversion of the NICER timestamps to absolute time in MJD can be done by 
		MJD_TT = MJDREFI+MJDREFF+(TIMEZERO+TIME)/ 86400


EXAMPLES
	1. Get the NICER TIME origion. 
		$ %prog 2014-01-01T00:00:00 -f isot -s utc

	2. Convert Calender format to the NICER time (126266402.000).
		$ %prog 2018-01-01T10:00:00 -f isot -s utc

	3. Convert NICER TIME to other time format (oposite of 2).
		$ %prog 126266402 -f met -s met 

	2. Convert MJD_UTC 58000 to other formats. 
		$ %prog 58000 -f mjd -s utc  

"""
parser = OptionParser(usage=usage)
#parser.add_option("-i","--input",dest="input",
 #       action="store",help="Input time value.",type="string")
parser.add_option("-f","--format",dest="format",default="met",
       action="store",help="Time format, any of (met, jd, mjd, isot)",type="string")
parser.add_option("-s","--scale",dest="scale",default="met",
       action="store",help="Time scale, any of (utc or tt)",type="string")
(options, args) = parser.parse_args()


input_value = args[0]

dump  = "----- Input Time Value and Formats -----\n"
dump += "intime: %s\n" % str(input_value)
dump += "format: %s\n" % options.format
dump += "scale : %s\n" % options.scale 

if options.format == "met":
	mission = float(input_value)
	mjd_tt  = NICER_MJDREFI+NICER_MJDREFF+(NICER_TIMEZERO+mission)/86400.0
	time_tt = Time(mjd_tt,format='mjd',scale='tt')
	time_utc = time_tt.utc
else:
	if input_value.isdigit():
		time = Time(float(input_value),format=options.format,scale=options.scale)
	else:
		time = Time(str(input_value),format=options.format,scale=options.scale)		
	time_tt = time.tt 
	time_tt.format = 'mjd'
	mjd_tt = time_tt 
	mission = (float(mjd_tt.mjd) - NICER_MJDREFI - NICER_MJDREFF) * 86400.0 - NICER_TIMEZERO
	time_utc = time.utc 
	time_utc.format = 'mjd'

fermi_time  = time_tt.gps - FERMI_MET_ORIGIN.gps
nustar_time = time_tt.gps - NUSTART_MET_ORIGIN.gps 
rxte_time   = time_tt.gps - RXTE_MET_ORIGIN.gps 
suzaku_time = time_tt.gps - SUZAKU_MET_ORIGIN.gps 
#swift_time  = time_tt.gps - SWIFT_MET_ORIGIN.gps 
xmm_time    = time_tt.gps - XMM_MET_ORIGIN.gps 
chandra_time = time_tt.gps - CHANDRA_MET_ORIGIN.gps 

dump += "----- Calendar Time Formats -----\n"
dump += "ISO8601_TT : %s (TT)\n" % time_tt.isot 
dump += "MJD_TT     :   %.8f (TT)\n" % time_tt.mjd
dump += " JD_TT     : %.8f (TT) \n" % time_tt.jd
dump += "ISO8601_UTC: %s (UTC)\n" % time_utc.isot 
dump += "MJD_UTC    :   %.8f (UTC) \n" % time_utc.mjd
dump += " JD_UTC    : %.8f (UTC) \n" % time_utc.jd
dump += "----- Mission-Specific Time Formats (Misson Elapsed Time, NET) -----\n"
dump += "Fermi seconds sicne 2001.0 UTC (decimal)     : %.8f\n" % fermi_time
dump += "NuSTAR seconds since 2010.0 UTC (decimal)    : %.8f\n" % nustar_time
dump += "RXTE seconds since 1994.0 UTC (decimal)      : %.8f\n" % rxte_time 
dump += "Suzaku seconds since 2000.0 UTC (decimal)    : %.8f\n" % suzaku_time
#dump += "Swift seconds since 2001.0 UTC (decimal): %.8f\n" % swift_time
dump += "XMM/Chandra seconds since 1998.0 TT (decimal): %.8f\n" % xmm_time
dump += "NICER seconds since 2014.0 UTC (decimal)     : %.8f" % mission 
print(dump)



