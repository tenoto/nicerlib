import os 
import sys 
import glob 
import numpy as np
import astropy.io.fits as pyfits

from pyheasoft import * 
from niconst import * 

NICER_DATA_SUBDIRECTORY_LIST = ['auxil','log','xti'] 

class NicerBackgroundModel():
	def __init__(self,obsid_path,outdir,exposure_threshold=60.0,
		rmffile="/Users/enoto/work/niresp/nicer_v1.02.rmf",
		arffile="/Users/enoto/work/niresp/ni_xrcall_onaxis_v1.02.arf"):
		self.obsid_path = obsid_path		
		self.outdir = outdir 
		self.exposure_threshold = exposure_threshold
		self.rmffile = rmffile
		self.arffile = arffile 
		sys.stdout.write('obsid_path : %s\n' % obsid_path)
		sys.stdout.write('outdir : %s\n' % outdir)
		sys.stdout.write('exposure_threshold : %.3f (s)\n' % exposure_threshold)		
		self.obsid = os.path.basename(self.obsid_path)

		self.bgddir = '%s/gtibgd' % self.outdir
		self.inputfile_bgd_param = '%s/ni%s_bgd_gti_param.txt' % (self.bgddir,self.obsid)
		self.inputfile_bgd_param_detail = self.inputfile_bgd_param.replace('.txt','_niprefilter2.txt')
		#print(self.inputfile_bgd_param)
		#print(self.inputfile_bgd_param_detail)

	def check_input_directory(self):
		sys.stdout.write('...Making an object NicerObservation...\n')
		if not os.path.exists(self.obsid_path):
			sys.stderr.write('Error: input directory does not exist.\n==> %s\n' % self.obsid_path)
			exit()
		else:
			sys.stdout.write('input direcotry is ok: %s\n' % self.obsid_path)
		for subdir in NICER_DATA_SUBDIRECTORY_LIST:
			if not os.path.exists('%s/%s' % (self.obsid_path,subdir)):
				sys.stderr.write('Error: input directory does not have the subdirectory [%s]\n' % subdir)				
				exit()
			else:
				sys.stdout.write('subdirectory %s exist.\n' % subdir)
		self.obsid = os.path.basename(self.obsid_path)

	def set_input_files(self):
		self.hk_dir      = '%s/xti/hk' % self.obsid_path 
		self.evtcl_dir   = '%s/xti/event_cl' % (self.obsid_path)
		self.evtuf_dir   = '%s/xti/event_uf' % (self.obsid_path)		

		self.attfile = '%s/auxil/ni%s.att.gz' % (self.obsid_path, self.obsid)
		if os.path.exists(self.attfile):
			sys.stdout.write('attfile: %s\n' % self.attfile)
		elif os.path.exists(self.attfile.replace('.gz','')):
			cmd = 'gzip %s' % self.attfile.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.attfile)
			exit()

		self.mkffile = '%s/auxil/ni%s.mkf.gz' % (self.obsid_path, self.obsid)
		if os.path.exists(self.mkffile):
			sys.stdout.write('mkffile: %s\n' % self.mkffile)
		elif os.path.exists(self.mkffile.replace('.gz','')):
			cmd = 'gzip %s' % self.mkffile.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.mkffile)
			exit()				

		self.orbfile = '%s/auxil/ni%s.orb.gz' % (self.obsid_path, self.obsid)	
		if os.path.exists(self.orbfile):
			sys.stdout.write('orbfile: %s\n' % self.orbfile)
		elif os.path.exists(self.orbfile.replace('.gz','')):
			cmd = 'gzip %s' % self.orbfile.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.orbfile)
			exit()			

		self.catfile = '%s/auxil/ni%s.cat.gz' % (self.obsid_path, self.obsid)	
		if os.path.exists(self.catfile):
			sys.stdout.write('catfile: %s\n' % self.catfile)
		elif os.path.exists(self.catfile.replace('.gz','')):
			cmd = 'gzip %s' % self.catfile.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.catfile)
			exit()		

		self.clevt_mpu7  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (self.obsid_path, self.obsid)
		if os.path.exists(self.clevt_mpu7):
			sys.stdout.write('clevt_mpu7: %s\n' % self.clevt_mpu7)
		elif os.path.exists(self.clevt_mpu7.replace('.gz','')):
			cmd = 'gzip %s' % self.clevt_mpu7.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.clevt_mpu7)
			exit()			

		self.ufaevt_mpu7 = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt.gz' % (self.obsid_path, self.obsid)		
		if os.path.exists(self.ufaevt_mpu7):
			sys.stdout.write('ufaevt_mpu7: %s\n' % self.ufaevt_mpu7)
		elif os.path.exists(self.ufaevt_mpu7.replace('.gz','')):
			cmd = 'gzip %s' % self.ufaevt_mpu7.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.ufaevt_mpu7)
			exit()				

		self.ufevt_mpu0 = '%s/xti/event_uf/ni%s_0mpu0_uf.evt.gz' % (self.obsid_path, self.obsid)		
		if os.path.exists(self.ufevt_mpu0):
			sys.stdout.write('ufevt_mpu0: %s\n' % self.ufevt_mpu0)
		elif os.path.exists(self.ufevt_mpu0.replace('.gz','')):
			cmd = 'gzip %s' % self.ufevt_mpu0.replace('.gz','')
			print(cmd);os.system(cmd)
		else:
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.ufevt_mpu0)
			exit()			
		
	def show_input_files(self):
		sys.stdout.write('---------------------------------\n')
		sys.stdout.write('obsid_path = %s\n' % self.obsid_path)
		sys.stdout.write('obsid = %s\n' % self.obsid)
		sys.stdout.write('hk directory = %s\n' % self.hk_dir)
		sys.stdout.write('event clean directory = %s\n' % self.evtcl_dir)
		sys.stdout.write('event uf directory = %s\n' % self.evtuf_dir)
		sys.stdout.write('attitude file = %s\n' % self.attfile)
		sys.stdout.write('mkf file = %s\n' % self.mkffile)						
		sys.stdout.write('orbit file = %s\n' % self.orbfile)								
		sys.stdout.write('clean event mpu7 file = %s\n' % self.clevt_mpu7)										
		sys.stdout.write('ufa event mpu7 file = %s\n' % self.ufaevt_mpu7)												
		sys.stdout.write('uf event mpu0 file = %s\n' % self.ufevt_mpu0)														
		sys.stdout.write('---------------------------------\n')

	def set_header_keywords(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		hdu = pyfits.open(self.ufevt_mpu0)
		self.date_obs = hdu['EVENTS'].header['DATE-OBS']
		self.date_end = hdu['EVENTS'].header['DATE-END']		
		self.tstart   = hdu['EVENTS'].header['TSTART']				
		self.tstop    = hdu['EVENTS'].header['TSTOP']						
		self.exposure = hdu['EVENTS'].header['EXPOSURE']		
		self.ontime   = hdu['EVENTS'].header['ONTIME']
		self.nevents  = hdu['EVENTS'].header['NAXIS2']	
		self.ra_nom   = hdu['EVENTS'].header['RA_NOM']
		self.dec_nom  = hdu['EVENTS'].header['DEC_NOM']
		self.ra_obj   = hdu['EVENTS'].header['RA_OBJ']
		self.dec_obj  = hdu['EVENTS'].header['DEC_OBJ']
		self.object   = hdu['EVENTS'].header['OBJECT']
		self.rate   = float(self.nevents)/float(self.exposure)		
		self.ufmpu0_total_gti_exposure = get_total_gti_exposure(self.ufevt_mpu0)

	def show_header_keywords(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		sys.stdout.write('DATE-OBS : %s\n' % self.date_obs)
		sys.stdout.write('DATE-END : %s\n' % self.date_end)
		sys.stdout.write('TSTART   : %s\n' % self.tstart)
		sys.stdout.write('TSTOP    : %s\n' % self.tstop)
		sys.stdout.write('EXPOSURE : %s\n' % self.exposure)
		sys.stdout.write('ONTIME   : %s\n' % self.ontime)
		sys.stdout.write('NEVENTS  : %s\n' % self.nevents)
		sys.stdout.write('RA_NOM   : %s\n' % self.ra_nom)
		sys.stdout.write('DECX_NOM : %s\n' % self.dec_nom)
		sys.stdout.write('RA_OBJ   : %s\n' % self.ra_obj)
		sys.stdout.write('DEC_OBJ  : %s\n' % self.dec_obj)		
		sys.stdout.write('OBJECT   : %s\n' % self.object)		
		sys.stdout.write('TOTAL_GTI_EXPOSURE : %.3f (s)\n' % self.ufmpu0_total_gti_exposure)

	def make_output_directory(self,outdir,flag_recreate=True):
		sys.stdout.write('=== %s (ObsID=%s) ===\n' % (sys._getframe().f_code.co_name,self.obsid))

		self.outdir = outdir
		if flag_recreate:
			cmd = 'rm -rf %s;mkdir -p %s' % (self.outdir,self.outdir)
		else:
			cmd = 'mkdir -p %s' % (self.outdir,self.outdir)			
		print(cmd);os.system(cmd)

	def make_bgdmodel_parameter_events(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		basename = os.path.basename(self.ufaevt_mpu7).replace('.gz','').replace('_ufa.evt','')

		self.xrayevt = '%s/%s_bx1x000.evt' % (self.outdir,basename)
		cmd  = 'rm -f %s;' % self.xrayevt
		cmd += 'fselect %s %s <<EOF\n' % (self.ufaevt_mpu7,self.xrayevt)
		cmd += 'EVENT_FLAGS==bx1x000\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)		

		self.IBGevt = '%s/%s_IBG.evt' % (self.outdir,basename)
		cmd  = 'rm -f %s;' % self.IBGevt
		cmd  = 'fselect %s %s <<EOF\n' % (self.xrayevt,self.IBGevt)
		cmd += '(PI >= 1500)&&(PI <= 1700)\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

		self.HREJevt = '%s/%s_HREJ.evt' % (self.outdir,basename)
		cmd  = 'rm -f %s;' % self.HREJevt
		cmd  = 'fselect %s %s <<EOF\n' % (self.xrayevt,self.HREJevt)
		cmd += '(PI >= 300)&&(PI <= 1800)&&(PI_RATIO > 1.54)\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

		self.NZevt = '%s/%s_NZ.evt' % (self.outdir,basename)
		cmd  = 'rm -f %s;' % self.NZevt
		cmd  = 'fselect %s %s <<EOF\n' % (self.xrayevt,self.NZevt)
		cmd += '(PI >= 0)&&(PI <= 200)\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

	def extract_spectrum(self,evtfile,gtifile,outpha):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		#outpha = evtfile.replace('.evt','.pha')
		cmd  = 'fxselect_extract_spec.py '
		cmd += '--inputevtfits=%s ' % evtfile
		cmd += '--outputpha=%s ' % outpha
		cmd += '--gtifile=%s ' % gtifile 
		cmd += '--rmf=%s ' % self.rmffile
		cmd += '--arf=%s ' % self.arffile
		print(cmd);os.system(cmd)

		return outpha

	def filter_gti_exposure(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		basename = os.path.basename(self.clevt_mpu7).replace('.gz','').replace('_cl.evt','')

		self.gtifile_expth = '%s/%s_gtisel%d.txt' % (self.outdir,basename,self.exposure_threshold)
		hdu_evt = pyfits.open(self.clevt_mpu7)
		f = open(self.gtifile_expth,'w')
		for line in hdu_evt['GTI'].data:
			exposure = float(line['STOP'])-float(line['START'])
			if exposure > self.exposure_threshold:
				dump = '%.6f,%.6f\n' % (line['START'],line['STOP'])
				f.write(dump)
		f.close()	

		outpha = self.gtifile_expth.replace('.txt','.pha')
		self.extract_spectrum(self.clevt_mpu7,self.gtifile_expth,outpha)
		#fxselect_extract_spec_gti(clevt,gtifile_expth,evtfile_expth)	

	def make_gti_bgd_parameters(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd = 'rm -rf %s; mkdir -p %s;' % (self.bgddir,self.bgddir)
		print(cmd);os.system(cmd)		

		mkfhdu = pyfits.open(self.mkffile)
		time = mkfhdu['PREFILTER'].data['TIME']
		gtinum = 1 
		f = open(self.inputfile_bgd_param,'w')
		f_detail = open(self.inputfile_bgd_param_detail,'w')
		for line in open(self.gtifile_expth):
			tmppha = '%s/ni%s_gti%04d.pha' % (self.bgddir,self.obsid,gtinum)
			#tmpevt = '%s/%s_gti%04d.evt' % (self.bgddir,outbase,gtinum)	
			tmpgti = tmppha.replace('.pha','.gti')
			f_tmp = open(tmpgti,'w')
			f_tmp.write(line)
			f_tmp.close()
			self.extract_spectrum(self.clevt_mpu7,tmpgti,tmppha)

			tstart = float(line.strip().split(',')[0])
			tstop  = float(line.strip().split(',')[1])

			cmd  = 'rm -f tmp.evt; '
			cmd += "fxselect_filter_time.py "
			cmd += "-i %s -o tmp.evt " % (self.IBGevt)
			cmd += "-d %.7f -u %.7f " % (tstart, tstop)
			print(cmd);os.system(cmd)
			tmphdu = pyfits.open("tmp.evt")
			nevents = len(tmphdu['EVENTS'].data)
			exposure = tmphdu['EVENTS'].header['EXPOSURE']
			IBG_rate = float(nevents)/float(exposure)			
			print(IBG_rate)

			cmd  = 'rm -f tmp.evt; '
			cmd += "fxselect_filter_time.py "
			cmd += "-i %s -o tmp.evt " % (self.HREJevt)
			cmd += "-d %.7f -u %.7f " % (tstart, tstop)
			print(cmd);os.system(cmd)
			tmphdu = pyfits.open("tmp.evt")
			nevents = len(tmphdu['EVENTS'].data)
			exposure = tmphdu['EVENTS'].header['EXPOSURE']
			HREJ_rate = float(nevents)/float(exposure)			
			print(HREJ_rate)

			cmd  = 'rm -f tmp.evt; '
			cmd += "fxselect_filter_time.py "
			cmd += "-i %s -o tmp.evt " % (self.NZevt)
			cmd += "-d %.7f -u %.7f " % (tstart, tstop)
			print(cmd);os.system(cmd)
			tmphdu = pyfits.open("tmp.evt")
			nevents = len(tmphdu['EVENTS'].data)
			exposure = tmphdu['EVENTS'].header['EXPOSURE']
			NZ_rate = float(nevents)/float(exposure)			
			print(NZ_rate)

			cmd  = 'rm -f tmp.evt; '
			print(cmd);os.system(cmd)

			gti_flag_array = np.where((tstart<=time)&(time<tstop),True,False)
			FPM_XRAY_PI_1500_1700_mean = mkfhdu['PREFILTER'].data['FPM_XRAY_PI_1500_1700'][gti_flag_array].mean()
			FPM_RATIO_REJ_COUNT_mean = mkfhdu['PREFILTER'].data['FPM_RATIO_REJ_COUNT'][gti_flag_array].mean()
			FPM_XRAY_PI_0000_0025_mean = mkfhdu['PREFILTER'].data['FPM_XRAY_PI_0000_0025'][gti_flag_array].mean()	
			FPM_NOISE25_COUNT_mean = mkfhdu['PREFILTER'].data['FPM_NOISE25_COUNT'][gti_flag_array].mean()			
			NUM_FPM_ON_mean = mkfhdu['PREFILTER'].data['NUM_FPM_ON'][gti_flag_array].mean()			
			dump  = "GTI %d (%.6f - %.6f)\n" % (gtinum,tstart,tstop)
			dump += "  Exposure = %d (s)\n" % len(time[gti_flag_array])
			dump += "  FPM_XRAY_PI_1500_1700_mean = %.6f (c/s)\n" % FPM_XRAY_PI_1500_1700_mean
			dump += "  FPM_RATIO_REJ_COUNT_mean = %.6f (c/s)\n" % FPM_RATIO_REJ_COUNT_mean
			dump += "  FPM_XRAY_PI_0000_0025_mean = %.6f (c/s)\n" % FPM_XRAY_PI_0000_0025_mean
			dump += "  FPM_NOISE25_COUNT_mean = %.6f (c/s)\n" % FPM_NOISE25_COUNT_mean
			dump += "  NUM_FPM_ON_mean = %d\n" % int(NUM_FPM_ON_mean)
			print(dump)
			f_detail.write(dump)

			dump = "%s %d %.6f %.6f %.6f\n" % (os.path.basename(tmppha),int(NUM_FPM_ON_mean),
				IBG_rate,HREJ_rate,NZ_rate)	
			print(dump)
			f.write(dump)
			gtinum += 1 

		f.close()
		f_detail.close()

	def make_gti_bgd_spectra(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		orgdir = os.getcwd()
		os.chdir(self.bgddir)
		modeldir = '%s/BGMod_3C50' % os.getenv('NICER_BGDMODEL_PATH')
		syslink_files = ['bg_group_3C50.table','day_noise_3C50.table','bg_group_3C50_ngt_*.pha','day_noise_nz*.pha','make_BG3C50_spectrum.go']
		for file in syslink_files:
			cmd = 'ln -s %s/%s .' % (modeldir,file)
			print(cmd);os.system(cmd)

		cmd  = 'rm -f setup_BG3C50.go;'
		cmd += 'awk \'{print "./make_BG3C50_spectrum.go", $1, $2, $3, $4, $5 " ;'
		cmd += 'cat calc_bg.go >> calc_BG3C50.xcm"}\' %s > setup_BG3C50.go' % os.path.basename(self.inputfile_bgd_param)
		print(cmd);os.system(cmd)

		cmd  = 'source setup_BG3C50.go'
		print(cmd);os.system(cmd)

		cmd  = 'source calc_BG3C50.xcm'
		print(cmd);os.system(cmd)

		for srcpha in glob.glob('ni%s_gti*.pha' % self.obsid):
			hdu = pyfits.open(srcpha)
			exposure = float(hdu['SPECTRUM'].header['EXPOSURE'])
			bgdpha = 'bg_%s' % srcpha
			#for extnum in range(3):
			for extnum in [1]:
				cmd = 'fparkey %.7f %s+%d EXPOSURE' % (exposure,bgdpha,extnum)
				print(cmd);os.system(cmd)

		for file in syslink_files:
			cmd = 'rm -f %s' % (file)
			print(cmd);os.system(cmd)			
		tentative_files = ['calc_bg.go','find_bg3C50_lib_spec.go','find_group.go','find_nzgroup.go','group.dat','inp.norms','nzgroup.dat','temp*','test*']
		for file in tentative_files:
			cmd = 'rm -f %s' % (file)
			print(cmd);os.system(cmd)			
		os.chdir(orgdir)

	def merge_bgd_spectra(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		orgdir = os.getcwd()
		os.chdir(self.bgddir)

		exposure_list = []
		for bgdpha in glob.glob('bg_*.pha'):
			hdu = pyfits.open(bgdpha)
			exposure_list.append(hdu['SPECTRUM'].header['EXPOSURE'])
		total_exposure = sum(exposure_list)

		expr = ""
		i = 0
		for bgdpha in glob.glob('bg_*.pha'):
			expr += "%s * %.6f / %.6f " % (bgdpha,exposure_list[i],total_exposure)
			i += 1 
			if i < len(glob.glob('bg_*.pha')):
				expr += "+ "
		print(expr)

		#outpha = 'bg_ni%s_gtisel%d.pha' % (self.obsid,self.exposure_threshold)
		outpha = '../bg_ni%s_0mpu7_gtisel%d.pha' % (self.obsid,self.exposure_threshold)
		cmd  = 'mathpha '
		cmd += '"%s" ' % expr
		cmd += 'R outfil="%s" ' % outpha
		cmd += 'exposure=%.6f ' % total_exposure
		cmd += 'errmeth=gaussian properr=yes ncomments=0 areascal=NULL clobber=yes'
		print(cmd);os.system(cmd)

		os.chdir(orgdir)

	def run(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.check_input_directory()
		self.set_input_files()
		self.show_input_files()
		self.set_header_keywords()
		self.show_header_keywords()
		self.make_output_directory(self.outdir,flag_recreate=True)
		self.make_bgdmodel_parameter_events()
		self.filter_gti_exposure()
		self.make_gti_bgd_parameters()
		self.make_gti_bgd_spectra()
		self.merge_bgd_spectra()






