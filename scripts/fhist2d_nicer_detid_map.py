#!/usr/bin/env python

__name__    = 'fhist2d_nicer_detid_map'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Mar. 29'

import os
import pyfits
import numpy as np 
import matplotlib.pyplot as plt
from optparse import OptionParser
from datetime import datetime 

from matplotlib import rc
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['font.size'] = 18
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.xmargin'] = '0' #'.05'
plt.rcParams['axes.ymargin'] = '.2'

RAWX_MAX = 7
RAWY_MAX = 6

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file",type="string")
parser.add_option("-o","--outputpdf",dest="outputpdf",
	action="store",help="output pdf file",type="string")
parser.add_option("-e","--expression",dest="expression",default=None, 
	action="store",help="filter expression",type="string")
parser.add_option("-t","--title",dest="title",
	action="store",help="title",type="string")
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "input fits file is needed. %> -i inputfits"
	quit()
if options.outputpdf == None:	
	print "output pdf file is needed. %> -o outputpdf"
	quit()
if options.expression == None:
	print "expression is need. %> -e expression"
	quit()
print "inputfits : %s " % options.inputfits
print "outputpdf: %s " % options.outputpdf
print "expression: %s " % options.expression
print "title: %s " % options.title

if not os.path.exists(options.inputfits):
	print "input file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputpdf):
	print "output pdf file has already existed: %s " % options.outputpdf
	quit()	

basename  = os.path.splitext(os.path.basename(options.outputpdf))[0]
f_tmpsel = 'tmpsel_%s.fits' % basename

cmd = 'rm -f %s' % (f_tmpsel)
print(cmd);os.system(cmd)

cmd = 'fselect %s %s "%s"' % (options.inputfits,f_tmpsel,options.expression)
print(cmd);os.system(cmd)

hist2d = np.zeros([RAWY_MAX+1,RAWX_MAX+1])
print(hist2d)

hdu_sel = pyfits.open(f_tmpsel)
array_rawx = hdu_sel['EVENTS'].data["RAWX"]
array_rawy = hdu_sel['EVENTS'].data["RAWY"]
for rawx in range(RAWX_MAX+1):
	for rawy in range(RAWY_MAX+1):
		array_rawxy_sel = np.logical_and(array_rawx==rawx,array_rawy==rawy)
		numsel_rawxy = len(array_rawx[array_rawxy_sel])
		hist2d[rawy,rawx] = numsel_rawxy
print(hist2d)

plt.gcf().clear()
fig = plt.figure(figsize=(14,10))
ax  = fig.add_subplot(111)
plt.xlabel('RAWX')
plt.ylabel('RAWY')
plt.title('title')
img = ax.imshow(hist2d,interpolation="nearest",origin="upper")
for i, cas in enumerate(hist2d):
	for j, c in enumerate(cas):
		if c>0:
			plt.text(j-.2, i+.2, c, fontsize=14)
cb = fig.colorbar(img)
fig.savefig(options.outputpdf,dpi=200)

cmd = 'rm -f %s' % (f_tmpsel)
print(cmd);os.system(cmd)

