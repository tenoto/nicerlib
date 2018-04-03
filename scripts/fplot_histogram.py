#!/usr/bin/env python

__name__    = 'fplot_histogram'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 April 2'

PLOT_DEV = '/dev/null'
#PLOT_DEV = '/xw'

import os 
import sys 
#cmd  = 'fhisto %s %s'
#fhisto infile outfile column binsz

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-x", "--logx", action="store_true", dest="logx", default=False)
parser.add_option("-y", "--logy", action="store_true", dest="logy", default=False)
parser.add_option("-f","--filter",dest="filter",default=None,
        action="store",help="filter expression",type="string")
(options, args) = parser.parse_args()

if len(args) != 6:
	sys.stderr.write('%s infits outpdf colname nbin min max\n' % __name__)
	quit()

infits = args[0]
outpdf = args[1]
colname = args[2]
nbin = int(args[3])
minv = float(args[4])
maxv = float(args[5])
print("infits: %s" % infits)
print("outpdf: %s" % outpdf)
print("colname: %s" % colname)
print("nbin: %d" % nbin)
print("minv: %.3e" % minv)
print("maxv: %.3e" % maxv)

ftmp_histo = '%s.fht' % os.path.splitext(outpdf)[0]
cmd = 'rm -f %s' % ftmp_histo
print(cmd);os.system(cmd)

if options.filter != None:
	tmp_infits = 'tmp_%s.fits' % os.path.splitext(os.path.basename(infits))[0]
	cmd = 'rm -f %s' % tmp_infits
	print(cmd);os.system(cmd)	
	cmd = 'fselect %s %s expr="%s" ' % (infits, tmp_infits, options.filter)
	print(cmd);os.system(cmd)

if options.filter != None:
	cmd  = 'fhisto infile=%s outfile=%s ' % (tmp_infits,ftmp_histo)
else:	
	cmd  = 'fhisto infile=%s outfile=%s ' % (infits,ftmp_histo)
cmd += 'column=%s binsz=%.6f lowval=%.6f highval=%.6f ' % (colname,float(maxv-minv)/float(nbin),minv,maxv)
cmd += 'outcolx=%s outcoly=NUMBERS ' % (colname)
print(cmd);os.system(cmd)

ftmp_ps = '%s.ps' % os.path.splitext(outpdf)[0]
cmd  = 'fplot %s %s NUMBERS - %s @ <<EOF\n' % (ftmp_histo,colname,PLOT_DEV)
cmd += 'la t %s\n' % infits
cmd += 'lwid 5 \n'
cmd += 'time off\n'
if options.logx:
	cmd += 'log x on\n'
if options.logy:
	cmd += 'log y on\n'	
cmd += 'line step on \n'
cmd += 'hard %s/cps\n' % ftmp_ps
cmd += 'exit\n'
cmd += 'EOF\n'
print(cmd);os.system(cmd)

cmd  = 'ps2pdf %s;\n' % ftmp_ps
if os.path.dirname(outpdf) != '':
	cmd += 'mv %s %s' % (os.path.basename(outpdf),os.path.dirname(outpdf))	
print(cmd);os.system(cmd)

cmd = 'rm -f %s %s' % (ftmp_histo,ftmp_ps)
print(cmd);os.system(cmd)

if options.filter != None:
	cmd = 'rm -f %s' % tmp_infits
	print(cmd);os.system(cmd)
