#!/usr/bin/env python

__name__    = 'fplot_scatter'
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

if len(args) != 4:
	sys.stderr.write('%s infits outpdf colx coly\n' % __name__)
	quit()

infits = args[0]
outpdf = args[1]
colx   = args[2]
coly   = args[3]
print("infits: %s" % infits)
print("outpdf: %s" % outpdf)
print("colx: %s" % colx)
print("coly: %s" % coly)

if options.filter != None:
	tmp_infits = 'tmp_%s.fits' % os.path.splitext(os.path.basename(infits))[0]
	cmd = 'rm -f %s' % tmp_infits
	print(cmd);os.system(cmd)	
	cmd = 'fselect %s %s expr="%s" ' % (infits, tmp_infits, options.filter)
	print(cmd);os.system(cmd)

ftmp_ps = '%s.ps' % os.path.splitext(outpdf)[0]
if options.filter != None:
	cmd  = 'fplot %s %s %s - %s @ <<EOF\n' % (tmp_infits,colx,coly,PLOT_DEV)
else:
	cmd  = 'fplot %s %s %s - %s @ <<EOF\n' % (infits,colx,coly,PLOT_DEV)
cmd += 'la t %s\n' % infits
cmd += 'lwid 5 \n'
cmd += 'mark 17 on 2\n'
cmd += 'time off\n'
if options.logx:
	cmd += 'log x on\n'
if options.logy:
	cmd += 'log y on\n'	
cmd += 'hard %s/cps\n' % ftmp_ps
cmd += 'exit\n'
cmd += 'EOF\n'
print(cmd);os.system(cmd)

cmd  = 'ps2pdf %s;\n' % ftmp_ps
if os.path.dirname(outpdf) != '':
	cmd += 'mv %s %s' % (os.path.basename(outpdf),os.path.dirname(outpdf))	
print(cmd);os.system(cmd)

cmd = 'rm -f %s' % (ftmp_ps)
print(cmd);os.system(cmd)

if options.filter:
	cmd = 'rm -f %s' % tmp_infits
	print(cmd);os.system(cmd)	

