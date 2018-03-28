#!/usr/bin/env python 

import sys
import os 
import subprocess

def get_nimissiontime_mjdtt(mjd_tt):
	"""
	aetimecalc version 2007-05-14
	Written by Y.ISHISAKI (TMU)
	Built on ANL HEADAS converter 1.81 for ANL version 1.81
	any of prompt/mission/date/date_tt/mjd/mjd_tt/yday[date_tt] date
	date string in UTC, 'yyyy-mm-ddThh:mm:ss.sss'[0000-00-00T00:00:00] 2014-01-01T00:00:00

	aetimecalc: *** show parameter ***

	       INPUT   'date'
	    LEAPFILE   '/usr/local/soft/heasoft/caldb/data/gen/bcf/leapsec_010117.fits' (CALDB;/usr/local/soft/heasoft/heasoft-6.22/x86_64-apple-darwin16.7.0/refdata/leapsec.fits)
	       DATE0   '2000-01-01T00:00:00.000' (51544.000 in MJD-UTC)
	       YDAY0   '2005-07-10T00:00:00.000' (53561.000 in MJD-UTC)
	     MJDREFI   51544
    	 MJDREFF   0.00074287037037037

	Mission Time = 441849603.000000 (s)
	Date in UTC  = 2014-01-01T00:00:00
	MJD  in UTC  = 56658.000000000 (dy)
	Date in TT   = 2014-01-01T00:01:07.183999
	MJD  in TT   = 56658.000777593 (dy)
	Y+3097.000 (dy)
	"""
	szk_mission_time0 = 441849603.000000

	cmd  = 'rm -f tmp_aetimecalc.log\n'
	cmd += 'aetimecalc input=mjd_tt mjd_tt=%s > tmp_aetimecalc.log\n' % mjd_tt
	os.system(cmd)

	cmd = 'grep "Mission Time" tmp_aetimecalc.log | awk \'{print $4}\''
	#szk_mission_time = float(commands.getoutput('grep "Mission Time" tmp_aetimecalc.log | awk \'{print $4}\''))
	szk_mission_time = float(subprocess.check_output(cmd,shell=True).split()[0])
	ni_mission_time = szk_mission_time - szk_mission_time0

	cmd  = 'rm -f tmp_aetimecalc.log\n'
	os.system(cmd)

	return ni_mission_time

if __name__ == "__main__":
	if len(sys.argv) != 2:
		sys.stderr.write('error: %s mjdtt\n' % sys.argv[0])
		quit()
	mjdtt = float(sys.argv[1])
	print("MJD_TT: %.7f" % mjdtt)
	print("NITIME: %.7f" % get_nimissiontime_mjdtt(mjdtt))

