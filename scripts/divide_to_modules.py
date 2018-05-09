#!/usr/bin/env python

import os 
import sys 
from nicerlib import niconst 

if len(sys.argv) != 2:
	sys.stderr.write('usage: %s infits' % sys.argv[0])
	quit()
inevt = sys.argv[1]
basename = inevt.replace('.gz','').replace('.evt','')

for mpu in range(0,niconst.NUM_OF_MPU):
	for fpm in range(0,niconst.NUM_OF_FPM):
		detid = "%02d" % (10*mpu+fpm)
		outevt = '%s_detid%s.evt' % (basename,detid)
		cmd = 'fselect %s %s expr="(DET_ID==%s)"' % (inevt,outevt,detid)
		print(cmd);os.system(cmd)
