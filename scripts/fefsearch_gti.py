#!/usr/bin/env python

__name__    = 'fefsearch_gti'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 March 20'

import os
import pyfits
from optparse import OptionParser
from datetime import datetime 

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file for efsearch",type="string")
parser.add_option("-o","--outputfits",dest="outputfits",
	action="store",help="output fits file",type="string")
parser.add_option("-p","--dper",dest="dper",
	action="store",help="Period (sec)",type="float")
parser.add_option("-d","--dpdot",dest="dpdot",default=0.0,
	action="store",help="Period derivative (sec/sec)",type="float")
parser.add_option("-n","--nphase",dest="nphase",
	action="store",help="Number  of  phases  in the folded light curve(s)",type="int")
parser.add_option("-r","--dres",dest="dres",
	action="store",help="period resolution in search (sec)",type="float")
parser.add_option("-s","--nper",dest="nper",
	action="store",help="number of period in search",type="int")
parser.add_option("-t","--threshold",dest="threshold",default=20,
	action="store",help="threshold cycle number",type="int")
parser.add_option("-f","--flagcleanup",dest="flagcleanup",default=True,
	action="store_false",help="flag clean up")
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "inputfits is needed. %> -i inputfits"
	quit()
if options.outputfits == None:
	print "outputfits is needed. %> -o outputfits"
	quit()	
if options.dper == None:
	print "dper is needed. %> -p dper"
	quit()
if options.dpdot == None:
	print "dpdot is needed. %> -d dpdot"
	quit()	
if options.nphase == None:
	print "nphase is needed. %> -n nphase"
	quit()
if options.dres == None:
	print "dres is needed. %> -r dres"
	quit()	
if options.nper == None:
	print "nper is needed. %> -s nper"
	quit()
print "inputfits : %s " % options.inputfits
print "outputfits: %s " % options.outputfits
print "dper: %.10f" % options.dper
print "dpdot: %.10e (s/s)" % options.dpdot
print "nphase: %d " % options.nphase
print "dper: %.10f (s)" % options.dper
print "nper: %d " % options.nper

if not os.path.exists(options.inputfits):
	print "input fits file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputfits):
	print "output fits file has already existed: %s " % options.outputfits
	quit()

threshold_exposure = options.dper * float(options.threshold)
outdir = 'tmpout_%s' % os.path.splitext(options.outputfits)[0]
cmd = 'rm -fr %s; mkdir -p %s' % (outdir,outdir)
os.system(cmd)

first_flag = True

hdu = pyfits.open(options.inputfits)
used_gti_index = []
used_gti_exp   = []
for i in range(len(hdu['GTI'].data)):
	tstart, tstop = hdu['GTI'].data[i]
	gtiexp = tstop - tstart
	print(i,tstart,tstop,gtiexp)

	if gtiexp <= threshold_exposure:
		continue
	used_gti_index.append(i)
	used_gti_exp.append(gtiexp)

	basename = os.path.splitext(os.path.basename(options.inputfits))[0]
	fout_evt = '%s/%s_%d.evt' % (outdir,basename,i)
	cmd  = 'fselect %s %s ' % (options.inputfits,fout_evt)
	cmd += 'expr="(TIME>=%.7f).AND.(TIME<%.7f)"' % (tstart,tstop)
	print(cmd); os.system(cmd)

	fout_efs = fout_evt.replace('.evt','.efs')
	cmd  = 'efsearch %s window="-" ' % fout_evt
	cmd += 'epochfo=1 sepoch=INDEF dper=%.7f nphase=%d nbint=INDEF ' % (options.dper,options.nphase)
	cmd += 'dres=%.7f nper=%d ' % (options.dres,options.nper)
	cmd += 'outfile=%s plot=no ' % fout_efs
	print(cmd); os.system(cmd)

	cmd  = 'fcalc infile=%s outfile=%s ' % (fout_efs,fout_efs)
	cmd += 'clname=CHI2_%d expr=CHISQRD1 clobber=yes \n' % i
	cmd += 'fdelcol %s+1 CHISQRD1 proceed=yes confirm=no\n' % fout_efs
	cmd += 'fdelcol %s+1 ERROR1 proceed=yes confirm=no\n' % fout_efs	
	print(cmd); os.system(cmd)

	if first_flag:
		cmd  = 'cp %s %s \n' % (fout_efs,options.outputfits)
		cmd += 'fcalc infile=%s outfile=%s ' % (options.outputfits,options.outputfits)
		cmd += 'clname=CHISQRDSUM expr=CHI2_%d clobber=yes \n' % i
		cmd += 'fdelcol %s+1 CHI2_%d proceed=yes confirm=no\n' % (options.outputfits,i)
		cmd += 'fparkey CHISDSUM %s+1 EXTNAME add=yes\n' % options.outputfits
		cmd += 'fparkey CHISDSUM %s+1 HDUCLAS3 add=yes\n' % options.outputfits		
		cmd += 'fappend %s+1 %s\n' % (fout_efs,options.outputfits)		
		cmd += 'fcalc infile=%s+1 outfile=%s ' % (options.outputfits,options.outputfits)
		cmd += 'clname=CHISQRDAVE expr=0.0 clobber=yes \n' 
		print(cmd); os.system(cmd)
		first_flag = False 
	else:
		cmd   = 'faddcol %s[2] %s CHI2_%d \n' % (options.outputfits,fout_efs,i)
		print(cmd); os.system(cmd)		
		
	if options.flagcleanup:
		cmd = 'rm -f %s %s' % (fout_evt, fout_efs)
		print(cmd); os.system(cmd)		

print used_gti_index
print used_gti_exp 
CHISQRDAVE_max_value  = 0.0
CHISQRDAVE_max_period = 0.0
hdu = pyfits.open(options.outputfits)
for i in range(len(hdu['CHISDSUM'].data)):
	#print("*****%d******" % i)
	chi2sum = 0.0
	for j in used_gti_index:
		chi2sum += hdu['RESULTS'].data[i]['CHI2_%d' % j]
		#print(hdu['RESULTS'].data[i]['CHI2_%d' % j])
	chi2ave = chi2sum/float(len(used_gti_index))
	hdu['CHISDSUM'].data[i]['CHISQRDSUM'] = chi2sum
	hdu['CHISDSUM'].data[i]['CHISQRDAVE'] = chi2ave
	if chi2ave >= CHISQRDAVE_max_value:
		CHISQRDAVE_max_value = chi2ave
		CHISQRDAVE_max_period = hdu['CHISDSUM'].data[i]['PERIOD']
hdu.writeto(options.outputfits,clobber=True)

cmd  = 'fparkey %.7f %s+1 BESTPERD add=yes\n' % (CHISQRDAVE_max_period,options.outputfits)
cmd += 'fparkey %.7f %s+1 MAXCHISD add=yes\n' % (CHISQRDAVE_max_value,options.outputfits)
print(cmd); os.system(cmd)	

fout_ps = options.outputfits.replace('.efs','.ps')
cmd  = 'fplot %s ' % options.outputfits
cmd += 'PERIOD CHISQRDAVE - /xw @ <<EOF\n'
cmd += 'time off\n'
cmd += 'lwid 5\n'
cmd += 'la f Max period = %.3f (chi2=%.3f)\n' % (CHISQRDAVE_max_period,CHISQRDAVE_max_value)
cmd += 'plot\n'
cmd += 'hard %s/cps\n' % fout_ps
cmd += 'exit\n'
cmd += 'EOF\n'
print(cmd); os.system(cmd)	
