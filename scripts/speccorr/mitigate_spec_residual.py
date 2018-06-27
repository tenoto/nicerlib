#!/usr/bin/env python

__name__    = 'mitigate_spec_residual'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 June 27'

PLOT_DEV = '/dev/null'

import os 
import sys 
import astropy.io.fits as pyfits 
from optparse import OptionParser

CRAB_DATA_DIR = '%s/scripts/speccorr/data' % os.getenv('NICER_SOFT_PATH')
CRAB_CORR_PHA = '%s/ni1013010110_crab_residual.pha' % CRAB_DATA_DIR
rmffile = '%s/nicer_v1.02.rmf' % os.getenv('NICER_RESP_PATH')
arffile = '%s/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf' % os.getenv('NICER_RESP_PATH')

parser = OptionParser()
parser.add_option("-i","--input_pha",dest="input_pha",default=None,
        action="store",help="input pha spectrum",type="string")
(options, args) = parser.parse_args()

if options.input_pha == None:
	sys.stderr.write('--input_pha is needed.\n')
	exit()

tmp_input_pha = 'tmp_%s' % os.path.basename(options.input_pha)
cmd = 'ln -s %s %s' % (options.input_pha, tmp_input_pha)
print(cmd);os.system(cmd)

cmd = 'ln -s %s .' % CRAB_CORR_PHA
print(cmd);os.system(cmd)
tmp_crab_corr_pha = os.path.basename(CRAB_CORR_PHA)

corrpha = '%s_crabcorr.pha'  % (os.path.splitext(os.path.basename(options.input_pha))[0])
cmd = 'rm -f %s' % corrpha
print(cmd);os.system(cmd)

hdu = pyfits.open(options.input_pha)
exposure = hdu['SPECTRUM'].header['EXPOSURE']

cmd  = "mathpha properr=yes errmeth='POISS-3' <<EOF\n" 
cmd += "%s/%s\n" % (tmp_input_pha,tmp_crab_corr_pha)
cmd += "R\n"
cmd += "%s\n" % corrpha
cmd += "%s\n" % exposure
cmd += "%\n"
cmd += "0\n"
cmd += "EOF\n"
print(cmd);os.system(cmd)

for hdunum in range(2):
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=RESPFILE \n" % (rmffile, corrpha, hdunum)
	print(cmd);os.system(cmd)
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=ANCRFILE \n" % (arffile, corrpha, hdunum)
	print(cmd);os.system(cmd)

