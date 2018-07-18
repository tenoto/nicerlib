#!/usr/bin/env python

__name__    = 'fcalc_countrate_deltat'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 July 17'

PLOT_DEV = '/dev/null'

import os 
import sys 
import numpy as np 
import astropy.io.fits as pyfits 
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-b","--binsz",dest="binsz",default=0.1,
	action="store",help="Binning size (sec).",type="float")
parser.add_option("-l","--lowval",dest="lowval",default=0.0,
	action="store",help="Histogram lower edge.",type="float")
parser.add_option("-u","--upval",dest="highval",default=10.0,
	action="store",help="Histogram higher edge.",type="float")
parser.add_option("--fit_xmin",dest="fit_xmin",default=None,
	action="store",help="fit range xmin.",type="float")
parser.add_option("--fit_xmax",dest="fit_xmax",default=None,
	action="store",help="fit range xmax.",type="float")
parser.add_option("--plot_xmin",dest="plot_xmin",default=0.0,
	action="store",help="Plot xmin.",type="float")
parser.add_option("--plot_xmax",dest="plot_xmax",default=10.0,
	action="store",help="Plot xmax.",type="float")
parser.add_option("--plot_ymin",dest="plot_ymin",default=0.1,
	action="store",help="Plot ymin.",type="float")
parser.add_option("--plot_ymax",dest="plot_ymax",default=1000.0,
	action="store",help="Plot ymax.",type="float")
(options, args) = parser.parse_args()

infits = args[0]
print("infits: %s" % infits)

binsz   = options.binsz 
lowval  = options.lowval
highval = options.highval
plot_xmin = options.plot_xmin
plot_xmax = options.plot_xmax
plot_ymin = options.plot_ymin
plot_ymax = options.plot_ymax

if not os.path.exists(infits):
	sys.stderr.write('error: file %s does not exist.\n' % infits)
	quit()


outpdf = '%s.pdf' % (os.path.splitext(os.path.basename(infits))[0])

title = '%s (bin size %f s)' % (os.path.basename(infits), binsz)

cmd = 'rm -f tmp1_delta_t.evt tmp2_delta_t.fht'
print(cmd);os.system(cmd)

cmd = 'fcalc %s tmp1_delta_t.evt DELTA_T "seqdiff(TIME)"' % infits
print(cmd);os.system(cmd)

hdu = pyfits.open(infits)
total_nevt = int(len(hdu['EVENTS'].data))
exposure = float(hdu['EVENTS'].header['EXPOSURE'])
lc_rate = float(total_nevt) / exposure
lc_rate_error = np.sqrt(float(total_nevt)) / exposure

cmd  = 'fhisto tmp1_delta_t.evt tmp2_delta_t.fht DELTA_T %.8f lowval=%.8f highval=%.8f ' % (binsz,lowval,highval)
cmd += 'outcolx="DELTA_T" outcoly="NEVENTS" extname="DELTA_T" '
print(cmd);os.system(cmd)

# ==================
# Plot and fitting 
# ==================

hdu = pyfits.open('tmp2_delta_t.fht')
print(hdu['DELTA_T'].data)

deltat  = hdu['DELTA_T'].data['DELTA_T']
nevents = hdu['DELTA_T'].data['NEVENTS']
nerror  = hdu['DELTA_T'].data['ERROR']

def fitFunc(x, a, b):
    return a*np.exp(-b*x) 

if options.fit_xmin != None and options.fit_xmax != None:
	flag = deltat>options.fit_xmin
	deltat_fit = deltat[flag]
	nevents_fit = nevents[flag]
else:
	deltat_fit = deltat
	nevents_fit = nevents

fitParams, fitCovariances = curve_fit(fitFunc, deltat_fit, nevents_fit)
print fitParams
print fitCovariances
sigma = [fitCovariances[0,0], fitCovariances[1,1]]
slope       = fitParams[1]
slope_error = sigma[1]
norm = fitParams[0]
norm_error = sigma[0]

slope2lc_ratio = lc_rate / slope 

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '0' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

# http://media.usm.maine.edu/~pauln/ScipyScriptRepo/CurveFitting.py
plt.figure(figsize=(6,5))
#plt.plot(deltat,nevents)
label1 = 'data (lc rate: %.3e+/-%.3e cps)' % (lc_rate,lc_rate_error)
plt.errorbar(deltat,nevents,fmt='bo',yerr=nerror,label=label1,zorder=1)
plt.title(title)
plt.yscale("log")
plt.xlabel("Delta_T (sec)")
plt.ylabel("Number of events / (%.5f ms)" % (binsz*1e+3))
plt.xlim(plot_xmin,plot_xmax)
plt.ylim(plot_ymin,plot_ymax)
label2 = 'exp fit (slope: %.3e+/-%.3e cps)' % (slope,slope_error)
plt.plot(deltat, fitFunc(deltat, fitParams[0], fitParams[1]),color='r',label=label2,zorder=2)
plt.legend()
plt.savefig(outpdf,bbox_inches=0,dpi=600)

cmd = 'rm -f tmp1_delta_t.evt tmp2_delta_t.fht'
print(cmd);os.system(cmd)

print('lc rate: %.4f +/- %.4f cps' % (lc_rate,lc_rate_error))
print('delta_t slope : %.4f +/- %.4f cps' % (slope,slope_error))
print('delta_t norm  : %.4e +/- %.4e' % (norm, norm_error))
print('slope2lc_ratio : %.4e ' % slope2lc_ratio)


