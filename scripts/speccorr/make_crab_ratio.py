#!/usr/bin/env python

__name__    = 'make_crab_ratio'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 June 27'

PLOT_DEV = '/dev/null'

import os 
import sys 
import astropy.io.fits as pyfits 
from optparse import OptionParser

CRAB_DATA_DIR = '%s/scripts/speccorr/data' % os.getenv('NICER_SOFT_PATH')
CRAB_SRC_PHA = '%s/ni1013010110_clmerge_screen_gtisel.pha' % CRAB_DATA_DIR
CRAB_BGD_PHA = '%s/ni1013010110_mitbgd_BGMod_3C50.pha' % CRAB_DATA_DIR
CRAB_SUB_PHA = '%s/ni1013010110_clmerge_screen_gtisel_bgdsub.pha' % CRAB_DATA_DIR
rmffile = '%s/nicer_v1.02.rmf' % os.getenv('NICER_RESP_PATH')
arffile = '%s/ni_xrcall_onaxis_v1.02_wo_detid14_34.arf' % os.getenv('NICER_RESP_PATH')

parser = OptionParser()
parser.add_option("-s","--srcpha",dest="srcpha",default=None,
        action="store",help="source pha spectrum",type="string")
parser.add_option("-b","--bgdpha",dest="bgdpha",default=None,
        action="store",help="background pha spectrum",type="string")
parser.add_option("-l","--emin",dest="emin",default=0.2,
        help="energy min",type="float")
parser.add_option("-u","--emax",dest="emax",default=10.0,
        help="energy max",type="float")
parser.add_option("-i","--min_significance",dest="min_significance",default=3,
        help="min significance",type="int")
parser.add_option("-n","--maxbin",dest="maxbin",default=30,
        help="max bin",type="int")
(options, args) = parser.parse_args()

if options.srcpha == None:
	sys.stderr.write('--srcpha is needed.\n')
	exit()
if options.bgdpha == None:
	sys.stderr.write('--bgdpha is needed.\n')
	exit()	

tmp_srcpha = 'tmp_%s' % os.path.basename(options.srcpha)
tmp_bgdpha = 'tmp_%s' % os.path.basename(options.bgdpha)
cmd = 'ln -s %s %s' % (options.srcpha, tmp_srcpha)
print(cmd);os.system(cmd)
cmd = 'ln -s %s %s' % (options.bgdpha, tmp_bgdpha)
print(cmd);os.system(cmd)

subpha = '%s_bgdsub.pha'  % (os.path.splitext(os.path.basename(options.srcpha))[0])
cmd = 'rm -f %s' % subpha
print(cmd);os.system(cmd)

hdu = pyfits.open(options.srcpha)
exposure = hdu['SPECTRUM'].header['EXPOSURE']

cmd  = "mathpha properr=yes errmeth='POISS-3' <<EOF\n" 
cmd += "%s-%s\n" % (tmp_srcpha,tmp_bgdpha)
cmd += "R\n"
cmd += "%s\n" % subpha
cmd += "%s\n" % exposure
cmd += "%\n"
cmd += "0\n"
cmd += "EOF\n"
print(cmd);os.system(cmd)

for hdunum in range(2):
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=RESPFILE \n" % (rmffile, subpha, hdunum)
	print(cmd);os.system(cmd)
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=ANCRFILE \n" % (arffile, subpha, hdunum)
	print(cmd);os.system(cmd)

cmd = 'ln -s %s .' % CRAB_SUB_PHA
print(cmd);os.system(cmd)
tmp_crabsubpha = os.path.basename(CRAB_SUB_PHA)

crabratio_pha = '%s_crabratio.pha' % (os.path.splitext(os.path.basename(options.srcpha))[0])
cmd  = 'rm -f %s' % crabratio_pha
print(cmd);os.system(cmd)

cmd  = "mathpha properr=yes errmeth='POISS-3' <<EOF\n" 
cmd += "%s/%s\n" % (subpha,tmp_crabsubpha)
cmd += "R\n"
cmd += "%s\n" % crabratio_pha
cmd += "%s\n" % exposure
cmd += "%\n"
cmd += "0\n"
cmd += "EOF\n"
print(cmd);os.system(cmd)

for hdunum in range(2):
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=RESPFILE \n" % (rmffile, crabratio_pha, hdunum)
	print(cmd);os.system(cmd)
	cmd = "fparkey value=%s fitsfile=%s+%d keyword=ANCRFILE \n" % (arffile, crabratio_pha, hdunum)
	print(cmd);os.system(cmd)

cmd = 'rm -f %s %s %s' % (tmp_srcpha,tmp_bgdpha,tmp_crabsubpha)
print(cmd);os.system(cmd)

crabratio_ps = crabratio_pha.replace('.pha','.ps')
crabratio_pdf = crabratio_pha.replace('.pha','.pdf')
cmd  = 'rm -f %s %s' % (crabratio_ps,crabratio_pdf)
print(cmd);os.system(cmd)

cmd  = 'xspec<<EOF\n'
cmd += 'data 1 %s\n' % crabratio_pha
cmd += 'setplot energy\n'
cmd += 'ignore **-%.2f %.2f-**\n' % (options.emin,options.emax)
cmd += 'setplot rebin %d %d\n' % (options.min_significance,options.maxbin)
cmd += 'iplot ld\n'
cmd += 'la t Crab Ratio\n'
cmd += 'la y Crab Ratio\n'
cmd += 'lwid 5\n'
cmd += 'lwid 5 on 1..100\n'
cmd += 'time off\n'
cmd += 'hard %s/cps\n' % crabratio_ps
cmd += 'exit\n'
print(cmd);os.system(cmd)

cmd = 'ps2pdf %s' % crabratio_ps
print(cmd);os.system(cmd)

cmd = 'rm -f %s' % crabratio_ps
print(cmd);os.system(cmd)














