import os
import sys 
import pyfits 

from niconst import * 

#####################################################
# General 
#####################################################

def ps2pdf(psfile):
	cmd  = 'ps2pdf %s\n' % psfile
	pdffile = '%s.pdf' % os.path.splitext(os.path.basename(psfile))[0]
	outdir = os.path.dirname(psfile)
	pdffile_path = '%s/%s' % (outdir,pdffile)
	cmd += 'mv %s %s\n' % (pdffile,pdffile_path)
	cmd += 'rm -f %s\n' % psfile
	os.system(cmd)
	sys.stdout.write(cmd)
	return pdffile_path	

#####################################################
# Fits file information 
#####################################################

def get_number_of_events(evt_fitsfile):
	if not os.path.exists(evt_fitsfile):
		sys.stderr.write('Error: file does not exist; %s' % evt_fitsfile)
	hdu = pyfits.open(evt_fitsfile)
	number_of_events = int(len(hdu['EVENTS'].data))
	return number_of_events

def get_number_of_gtis(evt_fitsfile):
	if not os.path.exists(evt_fitsfile):
		sys.stderr.write('Error: file does not exist; %s' % evt_fitsfile)
	hdu = pyfits.open(evt_fitsfile)
	number_of_gtis = int(len(hdu['GTI'].data))
	return number_of_gtis	

def get_total_gti_exposure(evtfile,flag_dump=True,extension_name='GTI'):
	hdu = pyfits.open(evtfile)

	gti_sum = 0.0
	if flag_dump:
		sys.stdout.write('GTI-num gti(s) sum(s) [start-stop]\n')
	for i in range(len(hdu[extension_name].data)):
		start = hdu[extension_name].data[i]['START']
		stop  = hdu[extension_name].data[i]['STOP']
		gti   = stop - start 
		gti_sum += gti
		dump = 'GTI-%d %.3f %.3f [%.3f -- %.3f]\n' % (i, gti, gti_sum, start,stop)
		if flag_dump:
			sys.stdout.write(dump)
	return gti_sum

def get_observation_date(evtfile,start_stop='start',extension_name='EVENTS'):
	hdu = pyfits.open(evtfile)

	if start_stop == 'start':
		return hdu[extension_name].header['DATE-OBS']
	elif start_stop == 'stop':
		return hdu[extension_name].header['DATE-END']	


def nicer_fselect_energy(infits,emin_keV,emax_keV):

	pi_min = int(KEV_TO_PI * emin_keV)
	pi_max = int(KEV_TO_PI * emax_keV)
	emin_str = str(emin_keV).replace('.','p')
	emax_str = str(emax_keV).replace('.','p')		

	if os.path.dirname(infits) == '':
		outdir = '.'
	else:
		outdir = os.path.dirname(infits)
	ext = os.path.splitext(os.path.basename(infits))[-1]	
	if ext == '.gz':
		outfits = '%s/%s_%s_%skeV.evt' % (
			outdir, 
			os.path.splitext(os.path.splitext(os.path.basename(infits))[0])[0],
			emin_str,emax_str)
	else:
		outfits = '%s/%s_%s_%skeV.evt' % (
			outdir, 
			os.path.splitext(os.path.basename(infits))[0],
			emin_str,emax_str)
	if os.path.exists(outfits):
		sys.stderr.write("event file is already created. skipped.")
	else:
		cmd  = 'fselect ' 
		cmd += '%s ' % infits
		cmd += '%s ' % outfits
		cmd += '"PI >= %d && PI < %d"' % (pi_min,pi_max)
		print(cmd); os.system(cmd)
	return outfits 

#####################################################
# Xselect wrapper 
#####################################################

def xselect_extract_curve(input_filelist,outlc,binsize,pi_min=None,pi_max=None):
	cmd = """
xselect <<EOF
xsel
read event 
./
%s 
yes
set binsize %d 
""" % (input_filelist, binsize)
	if pi_min!=None and pi_max!=None:
		cmd += """
filter pha_cutoff
%d
%d
""" % (pi_min,pi_max)
	cmd += """
show filter 	
show data
extract curve




save curve
%s
exit
no
EOF
"""  % outlc
	print(cmd);os.system(cmd)

	log = outlc.replace('.lc','_lc.log').replace('.flc','_flc.log')
	cmd = 'mv xselect.log %s\n' % log
	print(cmd);os.system(cmd)

def xselect_extract_spectrum(evtfile,outdir=None,rmffile=None,arffile=None,dict_keywords=None):
	sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

	if outdir == None:
		phafile = '%s.pha' % os.path.splitext(evtfile)[0]
		phafile_log = '%s_pha.log' % os.path.splitext(evtfile)[0]
	else:		
		phafile = '%s/%s.pha' % (outdir,os.path.basename(os.path.splitext(evtfile)[0]))
		phafile_log = '%s/%s_pha.log' % (outdir,os.path.basename(os.path.splitext(evtfile)[0]))
	"""		
	# This is the original, in order to include 
	# echo "NICER:trunc_spec_lo 0" >> $HEADAS/bin/xselect.mdb
	# xselect is used. 
	cmd = "extractor filename=%s " % (evtfile)
	cmd += 'eventsout=NONE imgfile=NONE phafile=%s ' % (phafile)
	cmd += 'fitsbinlc=NONE regionfile=NONE timefile=NONE xcolf=RAWX ycolf=RAWY tcol=TIME '
	cmd += 'ecol=PI xcolh=RAWX ycolh=RAWY gti=GTI >& %s ' % phafile_log
	"""			
	cmd  = 'xselect << EOF\n'
	cmd += 'no\n'
	cmd += 'xsel\n'
	cmd += 'read event %s ./\n'  % evtfile
	cmd += 'yes\n'
	cmd += 'set phaname PI\n'
	cmd += 'show data\n'
	cmd += 'show filter\n'
	cmd += 'extract spectrum\n'
	cmd += '\n'
	cmd += '\n'	
	cmd += '\n'
	cmd += '\n'		
	cmd += 'save spectrum\n'
	cmd += '%s\n' % phafile
	cmd += 'exit\n'
	cmd += 'no\n'
	cmd += 'EOF\n'
	cmd += 'mv xselect.log %s' % phafile_log
	print(cmd); os.system(cmd)				

	hdunum = 1 
	if rmffile != None:
		cmd = "fparkey value=%s fitsfile=%s+%d keyword=RESPFILE \n" % (rmffile, phafile, hdunum)
		print(cmd);os.system(cmd)			
	if arffile != None:		
		cmd = "fparkey value=%s fitsfile=%s+%d keyword=ANCRFILE \n" % (arffile, phafile, hdunum)
		print(cmd);os.system(cmd)			

	print(dict_keywords)
	hdu = pyfits.open(phafile)
	if dict_keywords != None:
		for key in dict_keywords:
			for hdunum in range(len(hdu)):
				cmd = 'fparkey value="%s" fitsfile=%s+%d keyword="%s" add=yes \n' % (dict_keywords[key], phafile, hdunum, key)
				print(cmd);os.system(cmd)				

	return phafile 	

#####################################################
# Plot fitsfile 
#####################################################	

#def plot_xspec_spectrum(phafile):
#	print("plot_xspec_spectrum")
	#outdir = '%s/spec' % self.outdir
	#self.phafile = extract_spectrum(self.merged_clevt,
	#		outdir=outdir, 
	#		rmffile=self.param['rmffile'],arffile=self.param['arffile'])

	"""
	self.xcmfile = '%s.xcm' % (os.path.splitext(self.phafile)[0])
	dump  = "data 1 %s\n" % self.phafile
	dump += "resp 1 %s/%s\n" % (os.getenv('NICER_RESP_PATH'),self.param['rmffile'])
	dump += "arf 1 %s/%s\n" % (os.getenv('NICER_RESP_PATH'),self.param['arffile'])
	dump += "setplot energy\n"
	dump += "ignore 1:**-0.4 15.0-**\n"
	dump += "setplot rebin %d %d\n" % (self.param['rebin_sigma'],self.param['rebin_num'])
	f = open(self.xcmfile,'w')
	f.write(dump)
	f.close()

	hdu = pyfits.open(self.phafile)
	self.pha_exposure = float(hdu[1].header['EXPOSURE'])
	self.pha_exposure_ks = self.pha_exposure * 1e-3
	self.pha_totcts   = int(hdu[1].header['TOTCTS'])
	self.pha_rate     = float(self.pha_totcts)/self.pha_exposure
	self.pha_date_obs = hdu[1].header['DATE-OBS']
	self.pha_date_end = hdu[1].header['DATE-END']
	self.pha_tstart = hdu[1].header['TSTART']
	self.pha_tstop = hdu[1].header['TSTOP']

	self.pha_title = '%s cnt / %.2f ks = %.2e cps' % (self.pha_totcts, self.pha_exposure_ks, self.pha_rate)
	#title = '%s (%d cnt / %.2f ks = %.2e cps)' % (self.obsid, self.pha_totcts, self.pha_exposure_ks, self.pha_rate)
	#subtitle = '%s -- %s' % (self.pha_date_obs, self.pha_date_end)
	outname = '%s_pha' % (os.path.splitext(self.phafile)[0])

	cmd  = "xspec <<EOF\n"
	cmd += "@%s\n" % self.xcmfile
	cmd += "ipl ld\n"
	cmd += "time off\n"
	cmd += "lab ot %s\n" % self.otitle 
	cmd += "lab t %s\n" % self.pha_title
	cmd += "lab f %s\n" % self.ftitle_time
	cmd += "lwid 5 \n"
	cmd += "lwid 5 on 1\n"
	cmd += "r y %s %s\n" % (self.param["flux_plot_ymin"],self.param["flux_plot_ymax"])
	cmd += "hard tmp.ps/cps\n"
	cmd += 'exit\n'
	cmd += 'exit\n'		
	cmd += '<<EOF\n'
	print(cmd); os.system(cmd)
	cmd = 'mv tmp.ps %s.ps\n'  % outname
	print(cmd); os.system(cmd)
	self.pdf_spectrum = ps2pdf('%s.ps' % outname)
	"""

def plot_rate_histogram(inflc,column="RATE",binsz=1,
		outcolx="OVERONRY_RATE",outcoly="NUMOFBIN"):

	outhst = '%s.fht' % os.path.splitext(inflc)[0]
	outps = '%s_fht.ps' % os.path.splitext(inflc)[0]	
	cmd  = 'rm -f %s %s' % (outhst,outps)
	print(cmd);os.system(cmd)

	cmd  = "fhisto infile=%s " % inflc
	cmd += "outfile=%s " % outhst
	cmd += "column=%s " % column
	cmd += "binsz=%d " % binsz
	cmd += "outcolx=%s outcoly=%s " % (outcolx,outcoly)
	print(cmd);os.system(cmd)

	cmd  = 'fplot %s ' % outhst
	cmd += 'OVERONRY_RATE NUMOFBIN - /xw @ <<EOF\n'
	cmd += 'log\n'
	cmd += 'r x 1 10000\n'
	cmd += 'line step on\n'
	cmd += 'lwid 5\n'
	cmd += 'lwid 5 on 1..100\n'
	cmd += 'time off\n'
	cmd += 'hard %s/cps\n' % outps
	cmd += 'exit\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)

	outpdf = ps2pdf(outps)
	return outpdf 

def plot_correlation(infile,colx,coly,outps,ext=1,pltcmd=""):

	cmd  = 'fplot %s[%d] ' % (infile,ext)
	cmd += '%s %s ' % (colx, coly)
	cmd += 'rows="-" device="/xw" pltcmd="@" <<EOF\n' 
	cmd += 'time off\n'
	cmd += 'lwid 5 \n'
	cmd += 'lwid 5 on 1..100\n'
	cmd += 'mark 1 on 2\n'
	cmd += '%s\n' % pltcmd 
	cmd += 'hard %s/cps\n' % outps
	cmd += 'quit\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)

	outpdf = ps2pdf(outps)
	return outpdf 

def plot_curve(inflc,colx="TIME",coly="RATE",colye="ERROR",ymin=1e-2,ymax=1e+3):
	sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

	psfile = '%s.ps' % os.path.splitext(inflc)[0]
	cmd  = 'fplot %s ' % inflc 
	cmd += 'xparm=%s yparm=%s[%s] ' % (colx,coly,colye)
	cmd += 'rows="-" device="/xw" pltcmd="@" <<EOF\n'
	cmd += 'log y on\n'
	cmd += 'r y %.3e %.3e\n' % (ymin,ymax)
	cmd += 'mark 17 on 2\n'
	cmd += 'lwid 5 \n'
	cmd += 'lwid 5 on 1..100\n'	
	cmd += 'time off\n'
	cmd += 'hard %s/cps\n' % psfile 
	cmd += 'exit\n'
	cmd += 'EOF\n'
	print(cmd);os.system(cmd)
	ps2pdf(psfile)

def plot_xspec_spectrum(inpha,rmf=None,arf=None,sigma=10,maxbin=50,
		emin=0.2,emax=15.0,ymin=0.01,ymax=1e+4,title="",otitle="",ftitle=""):
	sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

	if not os.path.exists(inpha):
		sys.stderr.write('file %s does not exist.' % inpha)
		quit()

	outdir   = os.path.dirname(inpha)
	basename = os.path.splitext(os.path.basename(inpha))[0]
	fxcm_read = '%s/%s_pha_read.xcm' % (outdir,basename)
	f = open(fxcm_read,'w')	
	dump  = "data 1 %s\n" % inpha
	if rmf != None:
		data += "resp %s\n" % rmf
	if arf != None:
		data += "arf %s \n" % arf 
	dump += "setplot energy\n"
	dump += "setplot rebin %d %d\n" % (sigma,maxbin)
	f.write(dump)
	f.close()

	fpco = '%s/%s_pha_read.pco' % (outdir,basename)
	f = open(fpco,'w')
	dump  = "r x %.1f %.1f\n" % (emin,emax)
	dump += "r y %.1e %.1e\n" % (ymin,ymax)
	dump += "lwid 5\n"
	dump += "lwid 5 on 1..100\n"
	dump += "time off\n"
	dump += "la ot %s\n" % otitle	
	dump += "la t %s\n" % title
	dump += "la f %s\n" % ftitle 
	f.write(dump)
	f.close()

	fps = '%s/%s_pha_read.ps' % (outdir,basename)
	cmd  = "xspec <<EOF\n"
	cmd += "@%s\n" % fxcm_read
	cmd += "ipl ld \n"
	cmd += "@%s\n" % fpco
	cmd += "hard %s/cps\n" % fps
	cmd += "quit\n"
	cmd += "EOF"
	os.system(cmd)

	return ps2pdf(fps)

#####################################################
# Timing 
#####################################################	

def get_nimissiontime_datett(date_tt):
	"""
	(nicer) [enoto@vicuna:ver0.06]$ aetimecalc
	aetimecalc version 2007-05-14
	Written by Y.ISHISAKI (TMU)
	Built on ANL HEADAS converter 1.81 for ANL version 1.81
	any of prompt/mission/date/date_tt/mjd/mjd_tt/yday[date_tt] 
	date string in TT, 'yyyy-mm-ddThh:mm:ss.sss'[2014-01-01T00:01:07.183989] 2014-01-01T00:00:00.0

	aetimecalc: *** show parameter ***

	       INPUT   'date_tt'
	    LEAPFILE   '/usr/local/soft/heasoft/heasoft-6.21/x86_64-apple-darwin16.6.0/refdata/leapsec.fits' (CALDB;/usr/local/soft/heasoft/heasoft-6.21/x86_64-apple-darwin16.6.0/refdata/leapsec.fits)
	       DATE0   '2000-01-01T00:00:00.000' (51544.000 in MJD-UTC)
	       YDAY0   '2005-07-10T00:00:00.000' (53561.000 in MJD-UTC)
	     MJDREFI   51544
	     MJDREFF   0.00074287037037037

	Mission Time = 441849535.816000 (s)
	Date in UTC  = 0000-00-00T00:00:00
	MJD  in UTC  = -678973.000000000 (dy)
	Date in TT   = 2014-01-01T00:00:00
	MJD  in TT   = 56658.000000000 (dy)
	Y-732534.000 (dy)
	"""
	szk_mission_time0 = 441849535.816000 

	cmd  = 'rm -f tmp_aetimecalc.log\n'
	cmd += 'aetimecalc input=date_tt Date_tt=%s > tmp_aetimecalc.log\n' % date_tt
	os.system(cmd)

	cmd = 'grep "Mission Time" tmp_aetimecalc.log | awk \'{print $4}\''
	#szk_mission_time = float(commands.getoutput('grep "Mission Time" tmp_aetimecalc.log | awk \'{print $4}\''))
	szk_mission_time = float(subprocess.check_output(cmd,shell=True).split()[0])
	ni_mission_time = szk_mission_time - szk_mission_time0

	cmd  = 'rm -f tmp_aetimecalc.log\n'
	os.system(cmd)

	return ni_mission_time

def barycentric_correction(evtfile,orbfile,ra,dec,ephem=""):
	sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

	outfile = "%s_bary.evt" % os.path.splitext(evtfile)[0]
	logfile = "%s.log" % os.path.splitext(outfile)[0]
	cmd  = 'barycorr ra=%.6f dec=%.6f ' % (ra,dec)
	cmd += 'infile=%s ' % evtfile
	cmd += 'outfile=%s ' % outfile
	cmd += 'orbitfiles=%s ' % orbfile
	if ephem == "":
		cmd += 'refframe=ICRS  ephem=JPLEPH.430 '
	else:
		cmd += '%s ' % ephem 
	cmd += ">& %s " % logfile 
	print(cmd); os.system(cmd)
	return outfile 		



"""
def characterise_lightcurve(flcfile, colname):
	cmd  = 'fstatistic %s %s - ' % (flcfile, colname)
	cmd += '| grep "maximum" | awk \'{print $7}\''
	#rate_max = float(commands.getoutput(cmd))
	rate_max = float(subprocess.check_output(cmd,shell=True).split()[0])

	cmd  = 'fstatistic %s %s - ' % (flcfile, colname)
	cmd += '| grep "minimum" | awk \'{print $7}\''
	#rate_min = float(commands.getoutput(cmd))
	rate_min = float(subprocess.check_output(cmd,shell=True).split()[0])

	cmd  = 'fstatistic %s %s - ' % (flcfile, colname)
	cmd += '| grep "mean" | awk \'{print $8}\''
	#rate_mean = float(commands.getoutput(cmd))
	rate_mean = float(subprocess.check_output(cmd,shell=True).split()[0])

	cmd  = 'fstatistic %s %s - ' % (flcfile, colname)
	cmd += '| grep "standard" | awk \'{print $9}\''
	#rate_std = float(commands.getoutput(cmd))
	rate_std = float(subprocess.check_output(cmd,shell=True).split()[0])

	return rate_max, rate_min, rate_mean, rate_std
"""		
