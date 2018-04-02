#!/usr/bin/env python

__name__    = 'fplot_hist2d_pulseprofile'
__author__  = 'Teruaki Enoto'
__version__ = '1.00'
__date__    = '2018 Apr. 2'

import os
import pyfits
import numpy as np 
import matplotlib.pyplot as plt
from optparse import OptionParser

from matplotlib.colors import LinearSegmentedColormap

def generate_cmap(colors):
	"""
	https://qiita.com/kenmatsu4/items/fe8a2f1c34c8d5676df8
	"""
	values = range(len(colors))
	vmax = np.ceil(np.max(values))
	color_list = []
	for v, c in zip(values, colors):
		color_list.append( ( v/ vmax, c) )
	return LinearSegmentedColormap.from_list('custom_cmap', color_list)

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

parser = OptionParser()
parser.add_option("-i","--inputfits",dest="inputfits",
	action="store",help="input fits file",type="string")
parser.add_option("-o","--outputpdf",dest="outputpdf",
	action="store",help="output pdf file",type="string")
parser.add_option("-d","--emin",dest="emin", default=0.1,
	action="store",help="Energy min (keV)",type="float")
parser.add_option("-u","--emax",dest="emax", default=15.0,
	action="store",help="Energy max (keV)",type="float")
parser.add_option("-n","--nbin_energy",dest="nbin_energy", default=40, 
	action="store",help="Number of bins for energy",type="int")
parser.add_option("-x","--phasemin",dest="phasemin", default=0.0,
	action="store",help="Phase min",type="float")
parser.add_option("-y","--phasemax",dest="phasemax", default=1.0,
	action="store",help="Phase max",type="float")
parser.add_option("-z","--nbin_phase",dest="nbin_phase", default=20, 
	action="store",help="Number of bins for phase",type="int")
(options, args) = parser.parse_args()

if options.inputfits == None:
	print "input fits file is needed. %> -i inputfits"
	quit()
if options.outputpdf == None:	
	print "output pdf file is needed. %> -o outputpdf"
	quit()
print "inputfits : %s " % options.inputfits
print "outputpdf: %s " % options.outputpdf

if not os.path.exists(options.inputfits):
	print "input file does not exists: %s" % options.inputfits
	quit()
if os.path.exists(options.outputpdf):
	print "output pdf file has already existed: %s " % options.outputpdf
	quit()	


outfits = '%s_2dhist.fits' % os.path.splitext(options.outputpdf)[0]
f_tmp1  = '%s_2dhist_tmp1.fits' % os.path.splitext(options.outputpdf)[0]
cmd = 'rm -f %s %s %s' % (outfits,f_tmp1,options.outputpdf)
print(cmd);os.system(cmd)

cmd = 'fcalc infile=%s outfile=%s clname=LOGENERGY expr="log10(PI/100.0)" ' % (options.inputfits,f_tmp1)
print(cmd);os.system(cmd)

xbinsz = float(options.phasemax - options.phasemin)/float(options.nbin_phase)
ybinsz = (np.log10(options.emax)-np.log10(options.emin))/float(options.nbin_energy)
cmd  = 'f2dhisto infile=%s outfil=%s ' % (f_tmp1,outfits)
cmd += 'xbinsz=%.3f ybinsz=%.3f ' % (xbinsz,ybinsz)
cmd += 'xcol=PHASE ycol=LOGENERGY '
cmd += 'xrange="%.3f,%.3f" ' % (options.phasemin,options.phasemax)
cmd += 'yrange="%.3f,%.3f" ' % (np.log10(options.emin),np.log10(options.emax))
print(cmd);os.system(cmd)

hdu = pyfits.open(outfits)
hist2d = hdu[0].data

hist2d_subratio = []
for phase_at_fixed_energy in hist2d:
	#print(phase_at_fixed_energy)
	sum_value = np.sum(phase_at_fixed_energy)
	if sum_value == 0:
		hist2d_subratio.append(list([0.0]*options.nbin_phase*2))
	else:
		ave = sum_value/float(options.nbin_phase)
		#print(sum_value,ave)
		phase_at_fixed_energy = (phase_at_fixed_energy-ave)/ave*100.0
		double_dump = list(phase_at_fixed_energy) + list(phase_at_fixed_energy)
		hist2d_subratio.append(double_dump)
#print(hist2d_subratio)	

for i in range(options.nbin_energy+1):
	logval = np.log10(options.emin) + ybinsz * i 
	energy = 10**logval
	print(i,logval,energy)
expr = "ENERGY:10**(%.3f+%.5f*Ny)" % (np.log10(options.emin),ybinsz)
print(expr)
cmd  = 'fparkey "%s" %s[0] ENRGEXPR add=yes\n' % (expr,outfits)
print(cmd);os.system(cmd)

cmd = 'rm -f %s' % (f_tmp1)
print(cmd);os.system(cmd)

plt.gcf().clear()
fig = plt.figure(figsize=(10,8))
ax  = fig.add_subplot(111)
plt.xlabel('Nbin of phase [two pulse cycle, Nphase=%d/cycle]' % options.nbin_phase )
plt.ylabel('Nbin of Log10(Energy/keV) [Emin=%.1f Emax=%.1f Nbin=%d]' % (options.emin,options.emax,options.nbin_energy) )
expr = "ENERGY=10**(%.3f+%.5f*Ny) keV" % (np.log10(options.emin),ybinsz)
plt.title(expr)
plt.axes().set_aspect('equal')
#plt.title(options.title)
#cm = generate_cmap(['#87CEEB', '#2E8B57', '#F4A460'])
cm = generate_cmap(['#00008B', '#aaaaab', '#FFFFFF', '#F4D793', '#F4A460'])
#img = ax.imshow(hist2d,interpolation="nearest",origin="lower",cmap=cm)
img = ax.imshow(hist2d_subratio,interpolation="nearest",origin="lower",cmap=cm)
cb = fig.colorbar(img)
cb.set_label('flucutation (%)', labelpad=-40, y=1.05, rotation=0)
fig.savefig(options.outputpdf,dpi=200)









