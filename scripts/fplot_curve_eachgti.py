#!/usr/bin/env python

__name__    = 'fplot_curve_eachgti'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Mar. 29'

import os
import pyfits
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file",type="string")
parser.add_option("-n","--outputbase",dest="outputbase",
	action="store",help="output basename",type="string")
parser.add_option("-o","--outputdir",dest="outputdir",
	action="store",help="output directory",type="string")
parser.add_option("-t","--gtithreshold",dest="gtithreshold",	
	action="store",help="Threshold exposure (sec) for each GTI",type="float")
parser.add_option("-b","--tbin",dest="tbin",
	action="store",help="Time binning (sec)",type="float")
parser.add_option("-d","--emin",dest="emin",
	action="store",help="Energy threshold min (keV)",type="float")
parser.add_option("-u","--emax",dest="emax",
	action="store",help="Energy threshold max (keV)",type="float")
parser.add_option("-x","--title",dest="title", default="GTI light curve",
	action="store",help="title for a plot",type="string")
parser.add_option("-y","--plotymin",dest="plotymin", default=0.1,
	action="store",help="plot ymin",type="float")
parser.add_option("-z","--plotymax",dest="plotymax", default=1e+5,
	action="store",help="plot ymax",type="float")
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "input fits file is needed. %> -i inputfits"
	quit()
if options.outputbase == None:
	print "output basename is needed. %> -n outputbase"
	quit()	
if options.outputdir == None:
	print "output directory is needed. %> -o outputdir"
	quit()		
if options.gtithreshold == None:
	print "GTI threshold is needed. %> -o gtithreshold"
	quit()	
if options.tbin == None:
	print "time binning is needed. %> -b tbin"
	quit()		
if options.emin == None:
	print "Emin is needed. %> -u emin"
	quit()					
if options.emax == None:
	print "Emax is needed. %> -u emax"
	quit()						
print "inputfits : %s " % options.inputfits
print "outputbase: %s " % options.outputbase
print "outputdir: %s " % options.outputdir
print "gtithreshold: %.3f (sec) " % options.gtithreshold
print "tbin: %.3f (sec) " % options.tbin
print "emin: %.3f (keV) " % options.emin
print "emax: %.3f (keV) " % options.emax

if not os.path.exists(options.inputfits):
	print "input file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputdir):
	print "output directory has already existed: %s " % options.outputdir
	quit()	

cmd = 'rm -rf %s; mkdir -p %s' % (options.outputdir,options.outputdir)
print(cmd);os.system(cmd)

hdu = pyfits.open(options.inputfits)

gti_trash = 0.0
gti_used  = 0.0
gtinum = 1

for gti in hdu['GTI'].data:
	gti_exposure = gti['STOP'] - gti['START']	
	if gti_exposure > options.gtithreshold:
		print(gtinum,gti['START'],gti['STOP'],gti_exposure)
		subdir = '%s/gti%d' % (options.outputdir,gtinum)
		cmd = 'rm -rf %s; mkdir -p %s' % (subdir,subdir)
		print(cmd);os.system(cmd)

		outevt = '%s/%s_gti%d.evt' % (subdir,options.outputbase,gtinum)
		cmd  = 'fxselect_filter_time.py -i %s -o %s ' % (options.inputfits, outevt)
		cmd += '-d %.7f -u %.7f' % (gti['START'],gti['STOP'])
		print(cmd);os.system(cmd)		

		outflc = '%s/%s_gti%d.flc' % (subdir,options.outputbase,gtinum)
		cmd  = 'fxselect_extract_curve.py -i %s -o %s ' % (outevt, outflc)
		cmd += '-t %.3f -d %.2f -u %.2f'  % (options.tbin, options.emin, options.emax)
		print(cmd);os.system(cmd)		

		start_mjd_tt = hdu['GTI'].header['MJDREFI'] + hdu['GTI'].header['MJDREFF'] + (hdu['GTI'].header['TIMEZERO']+gti['START'])/86400.0
		stop_mjd_tt = hdu['GTI'].header['MJDREFI'] + hdu['GTI'].header['MJDREFF'] + (hdu['GTI'].header['TIMEZERO']+gti['STOP'])/86400.0		

		outps = '%s/%s_gti%d_flc.ps' % (subdir,options.outputbase,gtinum)
		cmd  = 'fplot %s TIME "RATE[ERROR]" - /xw @<<EOF\n' % outflc
		cmd += 'line step on \n'
		cmd += 'log y on\n'
		cmd += 'r y %.3e %.3e\n' % (options.plotymin,options.plotymax)
		cmd += 'lwid 5\n'
		cmd += 'time off\n'
		cmd += 'la t %s GTI-%d (%.2f-%.2f)\n' % (options.title,gtinum,gti['START'],gti['STOP'])
		cmd += 'la f %d-s bin, %.1f-%.1f keV, MJD_TT %.4f-%.4f\n' % (options.tbin, options.emin, options.emax, start_mjd_tt, stop_mjd_tt)
		cmd += 'hard %s/cps\n' % outps
		cmd += 'exit\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)		

		outpdf = '%s.pdf' % os.path.splitext(os.path.basename(outps))[0]
		cmd  = 'ps2pdf %s\n' % outps
		cmd += 'mv %s %s\n' % (outpdf,subdir)
		print(cmd);os.system(cmd)		

		gti_used += gti_exposure
	else:
		gti_trash += gti_exposure
	gtinum += 1 
print("gti_trash: %.02f" % gti_trash)		
print("gti_used: %.02f" % gti_used)		










