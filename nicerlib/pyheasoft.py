import os
import sys 
import pyfits 

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

def extract_spectrum(evtfile,outdir=None,rmffile=None,arffile=None):
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
	cmd += 'xsel\n'
	cmd += 'read event %s .\n'  % evtfile
	cmd += 'yes\n'
	cmd += 'set phaname PI\n'
	cmd += 'extract spectrum\n'
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

	return phafile 	

def plot_xspec_spectrum(phafile):
	print("hoge")
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
