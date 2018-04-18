import os 
import sys 
import yaml
import glob 
import pyfits

from pyheasoft import * 
from niconst import * 

NICER_DATA_SUBDIRECTORY_LIST = ['auxil','log','xti'] 
NICER_OUTPUT_SUBDIRS = ['proc/gti','timing']
NICER_YAML_KEYWORDS = ["title","outbase","overonly_lc_binsize","overonly_rate_thresholds","rmffile","arffile","ra","dec","ephem","pulse_search_overonly_threshold","pulse_search_emin_keV","pulse_search_emax_keV","lc_emin_keV","lc_emax_keV","lc_binsize"]

class NicerObservation():
	def __init__(self,obsid_path):
		self.obsid_path = obsid_path
		self.check_input_directory()
		self.set_input_files()
		self.show_input_files()
		self.set_header_keywords()
		self.show_header_keywords()

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
		if not os.path.exists(self.attfile):
			self.attfile = '%s/auxil/ni%s.att' % (self.obsid_path, self.obsid)
		elif not os.path.exists(self.attfile):
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.attfile)
			exit()

		self.mkffile = '%s/auxil/ni%s.mkf.gz' % (self.obsid_path, self.obsid)
		if not os.path.exists(self.mkffile):
			self.mkffile = '%s/auxil/ni%s.mkf' % (self.obsid_path, self.obsid)
		elif not os.path.exists(self.mkffile):
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.mkffile)		
			exit()

		self.orbfile = '%s/auxil/ni%s.orb.gz' % (self.obsid_path, self.obsid)	
		if not os.path.exists(self.orbfile):
			self.orbfile = '%s/auxil/ni%s.orb' % (self.obsid_path, self.obsid)	
		elif not os.path.exists(self.orbfile):
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.orbfile)		
			exit()

		self.catfile = '%s/auxil/ni%s.cat.gz' % (self.obsid_path, self.obsid)	
		if not os.path.exists(self.catfile):
			self.orbfile = '%s/auxil/ni%s.cat' % (self.obsid_path, self.obsid)	
		elif not os.path.exists(self.catfile):
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.catfile)		
			exit()			

		self.clevt_mpu7  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (self.obsid_path, self.obsid)
		if not os.path.exists(self.clevt_mpu7):
			self.clevt_mpu7  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt' % (self.obsid_path, self.obsid)
		elif not os.path.exists(self.clevt_mpu7):
			sys.stderr.write('Warning: file does not exist.\n==>%s' % self.clevt_mpu7)		
			self.clevt_mpu7 = None

		self.ufaevt_mpu7 = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt.gz' % (self.obsid_path, self.obsid)		
		if not os.path.exists(self.ufaevt_mpu7):
			self.ufaevt_mpu7 = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt' % (self.obsid_path, self.obsid)					
		elif not os.path.exists(self.clevt_mpu7):
			sys.stderr.write('Error: file does not exist.\n==>%s' % self.ufaevt_mpu7)		
			self.ufaevt_mpu7 = None

		self.ufevt_mpu0 = '%s/xti/event_uf/ni%s_0mpu0_uf.evt.gz' % (self.obsid_path, self.obsid)		
		if not os.path.exists(self.ufevt_mpu0):
			self.ufevt_mpu0 = '%s/xti/event_cl/ni%s_0mpu0_uf.evt' % (self.obsid_path, self.obsid)					
		elif not os.path.exists(self.ufevt_mpu0):
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

	def run_nicerl2(self,flag_skip=False):
		sys.stdout.write('=== %s (ObsID=%s) ===\n' % (sys._getframe().f_code.co_name,self.obsid))

		if flag_skip:
			print("skip of run_nicerl2")
		else:
			if os.path.exists(self.evtcl_dir):
				cmd = 'rm -f %s/*' % self.evtcl_dir
				print(cmd);os.system(cmd)

			cmd = 'nicerl2 %s' % self.obsid_path
			print(cmd);os.system(cmd)	

			cmd  = 'gzip %s\n' % (self.clevt_mpu7.replace('.gz',''))
			cmd += 'gzip %s\n' % (self.ufaevt_mpu7.replace('.gz',''))
			cmd += 'gzip %s\n' % (self.mkffile.replace('.gz',''))			
			cmd += 'gzip %s\n' % (self.catfile.replace('.gz',''))						
			print(cmd);os.system(cmd)

			self.clevt_mpu7  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (self.obsid_path, self.obsid)
			self.ufaevt_mpu7 = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt.gz' % (self.obsid_path, self.obsid)		

	def run_niprefilter2(self):
		sys.stdout.write('=== %s (ObsID=%s) ===\n' % (sys._getframe().f_code.co_name,self.obsid))

		niprefilter2_outfile = '%s_lv2.mkf' % os.path.splitext(self.mkffile.replace('.gz',''))[0]
		cmd  = 'rm -f %s.gz %s\n' % (niprefilter2_outfile,niprefilter2_outfile)
		print(cmd);os.system(cmd)

		cmd  = 'niprefilter2 '
		cmd += 'indir=%s ' % self.obsid_path
		cmd += 'infile=%s ' % self.mkffile
		cmd += 'outfile=%s  ' % niprefilter2_outfile
		print(cmd);os.system(cmd)

		cmd  = 'gzip %s\n' % niprefilter2_outfile
		print(cmd);os.system(cmd)

		self.mkffile = '%s.gz' % niprefilter2_outfile

	"""
	def extract_overonly_event(self,flag_clean=True):
		'Extract OVERHOOT Only count'
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		expr = '(EVENT_FLAGS==bxxx010)'
		fname_inputfiles = '%s/ni%s_0mpu7_overonly.lst' % (self.evtcl_dir,self.obsid)
		f = open(fname_inputfiles,'w')
		for mpunum in range(NUM_OF_MPU):
			uffile = '%s/ni%s_0mpu%d_uf.evt.gz[%s]' % (self.evtuf_dir,self.obsid,mpunum,expr)
			f.write('%s\n' % uffile)
		f.close()

		outevt = '%s/tmp_ni%s_0mpu7_overonly.evt' % (self.evtcl_dir,self.obsid)
		cmd = 'ftmerge @%s %s clobber=yes\n' % (fname_inputfiles,outevt)
		print(cmd);os.system(cmd)

		self.overonly_evt = '%s/ni%s_0mpu7_overonly.evt' % (self.evtcl_dir,self.obsid)
		cmd = 'niextract-events %s %s timefile=%s[GTI] clobber=yes\n' % (outevt,self.overonly_evt,self.clevt_mpu7)
		print(cmd);os.system(cmd)		

		if flag_clean:
			cmd = 'rm -f %s' % outevt
			print(cmd);os.system(cmd)

		cmd = 'ftcalc %s %s PI 0 clobber=yes' % (self.overonly_evt,self.overonly_evt)
		print(cmd);os.system(cmd)		

		for i in range(3):
			cmd  = 'fparkey %d %s+%d TLMIN10 add=yes\n' % (NICER_PI_MIN,self.overonly_evt,i)
			cmd += 'fparkey %d %s+%d TLMAX10 add=yes\n' % (NICER_PI_MAX,self.overonly_evt,i)			
			print(cmd);os.system(cmd)		

	def extract_energy_selected_xrayevent(self,emin,emax,flag_clean=True):
		'Extract energy-selected X-ray events'
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		expr = '(EVENT_FLAGS==bx1x000)'
		fname_inputfiles = '%s/ni%s_0mpu7_eselxray.lst' % (self.evtcl_dir,self.obsid)
		f = open(fname_inputfiles,'w')
		for mpunum in range(NUM_OF_MPU):
			uffile = '%s/ni%s_0mpu%d_uf.evt.gz[%s]' % (self.evtuf_dir,self.obsid,mpunum,expr)
			f.write('%s\n' % uffile)
		f.close()

		outevt = '%s/tmp_ni%s_0mpu7_eselxray.evt' % (self.evtcl_dir,self.obsid)
		cmd = 'ftmerge @%s %s clobber=yes\n' % (fname_inputfiles,outevt)
		print(cmd);os.system(cmd)

		self.eselxray_evt = '%s/ni%s_0mpu7_eselxray.evt' % (self.evtcl_dir,self.obsid)
		cmd = 'niextract-events %s %s timefile=%s[GTI] clobber=yes\n' % (outevt,self.eselxray_evt,self.clevt_mpu7)
		print(cmd);os.system(cmd)		

		if flag_clean:
			cmd = 'rm -f %s' % outevt
			print(cmd);os.system(cmd)

		cmd = 'ftcalc %s %s PI 0 clobber=yes' % (self.overonly_evt,self.overonly_evt)
		print(cmd);os.system(cmd)		

		for i in range(3):
			cmd  = 'fparkey %d %s+%d TLMIN10 add=yes\n' % (NICER_PI_MIN,self.overonly_evt,i)
			cmd += 'fparkey %d %s+%d TLMAX10 add=yes\n' % (NICER_PI_MAX,self.overonly_evt,i)			
			print(cmd);os.system(cmd)	
	"""

	def run_barycentric_correction(self,ra,dec,ephem=""):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.clevt_bary = barycentric_correction(self.clevt_mpu7,self.orbfile,ra=ra,dec=dec,ephem=ephem)


class NicerProcess():
	def __init__(self,args):
		if args.obsid_path[0] == '@':
			self.flag_multimode = True
			self.filename_obsid_path = args.obsid_path.split('@')[-1]
		else:
			self.flag_multimode = False
			self.obsid_path = args.obsid_path
		self.fparam = args.fparam
		self.outdir = args.outdir
		self.title  = args.title 
		self.outbase = args.outbase

	def show_input_parameters(self):
		sys.stdout.write("flag_multimode = %s\n" % self.flag_multimode)
		sys.stdout.write("fparam = %s\n" % self.fparam)
		sys.stdout.write("outdir = %s\n" % self.outdir)
		sys.stdout.write("title = %s\n" % self.title)
		sys.stdout.write("outbase = %s\n" % self.outbase) 

	def set_obsid_path_list(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.obsid_path_list = []
		if self.flag_multimode:
			sys.stdout.write('filename_obsid_path: %s\n' % self.filename_obsid_path)
			for line in open(self.filename_obsid_path):
				cols = line.split()
				if cols[0] == '#':
					continue
				if cols[0][-1] == '/':
					path = cols[0][:-1]
				else:
					path = cols[0]
				self.obsid_path_list.append(path)
		else:
			if self.obsid_path[-1] == '/':
				path = self.obsid_path[:-1]
			else:
				path = self.obsid_path		
			self.obsid_path_list.append(path)
		print(self.obsid_path_list)

	def show_obsid_path_list(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		print(self.obsid_path_list)

	def make_output_directory(self,flag_recreate=True):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		if (not os.path.exists(self.outdir)) or flag_recreate:
			cmd  = 'rm -rf %s; ' % self.outdir
			for subdir in NICER_OUTPUT_SUBDIRS:
				cmd += 'mkdir -p %s/%s;' % (self.outdir,subdir)
			print(cmd); os.system(cmd)
		elif flag_recreate==False and os.path.exists(self.outdir):
			sys.stdout.write('warning: directory has already existed.\n')

	def load_parameterfile(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		if not os.path.exists(self.fparam):
			sys.stderr.write('parameter yaml file does not exist. %s\n' % self.fparam)
			exit()
		self.param = yaml.load(open(self.fparam))
		sys.stdout.write('parameter file %s is loaded.\n' % self.fparam)

		if self.title == None:
			if "title" in self.param:
				self.title = self.param['title']
			else:
				self.title = "no title"
				self.param['title'] = "no title"
		else:
			self.param['title'] = self.title 

		if self.outbase == None:
			if "outbase" in self.param:
				self.outbase = self.param['outbase']
			else:
				self.outbase = "niout"
				self.param['outbase'] = "niout"
		else:
			self.param['outbase'] = self.outbase 

		self.check_yaml_keywords()

	def check_yaml_keywords(self):
		for key in NICER_YAML_KEYWORDS:
			if not self.param.has_key(key):
				sys.stdout.write('yaml file does not have a keyword %s.\n' % key)
				quit()

	def set_nicer_observations(self):
		self.niobs_list = []
		for obsid_path in self.obsid_path_list:
			self.niobs_list.append(NicerObservation(obsid_path))

		self.tstart = self.niobs_list[0].tstart
		self.tstop  = self.niobs_list[-1].tstop
		self.date_obs = self.niobs_list[0].date_obs
		self.date_end = self.niobs_list[-1].date_end
		self.ra_obj = self.niobs_list[0].ra_obj
		self.dec_obj = self.niobs_list[0].dec_obj		
		self.ra_nom = self.niobs_list[0].ra_nom
		self.dec_nom = self.niobs_list[0].dec_nom

	def run_nicerl2(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		for niobs in self.niobs_list:
			niobs.run_nicerl2()

	def run_niprefilter2(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		for niobs in self.niobs_list:
			niobs.run_niprefilter2()

	"""
	def extract_overonly_event(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		for niobs in self.niobs_list:
			niobs.extract_overonly_event()
	"""

	"""
	def set_title(self):
		if len(self.niobs_list) == 1:
			self.otitle = '%s (%s) process.v=%s' % (self.param['title'],self.niobs_list[0].obsid,__version__)
		else:
			self.otitle = '%s (%s...%s) [%d obsids] process.v=%s' % (self.param['title'],
				self.niobs_list[0].obsid,self.niobs_list[-1].obsid,len(self.niobs_list),__version__)
		self.ftitle_time = '%s - %s (%s - %s)' % (self.date_obs,self.date_end,self.tstart,self.tstop)
	"""

	def merge_mkffiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.mkflist = '%s/proc/%s_mgd_mkf.lst' % (self.outdir,self.param['outbase'])
		f = open(self.mkflist,'w')
		for niobs in self.niobs_list:
			sys.stdout.write('%s\n' % niobs.mkffile)
			f.write('%s\n' % niobs.mkffile)
		f.close()

		self.merged_mkffile = '%s/proc/%s_mgd.mkf' % (self.outdir,self.param['outbase'])
		cmd = 'ftmerge infile=@%s outfile=%s' % (self.mkflist,self.merged_mkffile)
		print(cmd); os.system(cmd)

	def merge_orbfiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.orblist = '%s/proc/%s_mgd_orb.lst' % (self.outdir,self.param['outbase'])
		f = open(self.orblist,'w')
		for niobs in self.niobs_list:
			sys.stdout.write('%s\n' % niobs.orbfile)
			f.write('%s\n' % niobs.orbfile)
		f.close()

		self.merged_orbfile = '%s/proc/%s_mgd.orb' % (self.outdir,self.param['outbase'])
		cmd = 'ftmerge infile=@%s outfile=%s' % (self.orblist,self.merged_orbfile)
		print(cmd); os.system(cmd)		

	"""
	def extract_overonly_curve(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.overonlylist = '%s/proc/%s_mgd_overonly.lst' % (self.outdir,self.param['outbase'])
		f = open(self.overonlylist,'w')		
		for niobs in self.niobs_list:
			sys.stdout.write('%s\n' % niobs.overonly_evt)
			f.write('%s\n' % niobs.overonly_evt)
		f.close()

		self.overonly_curve = '%s/proc/%s_mgd_overonly.flc' % (self.outdir,self.param['outbase'])
		xselect_extract_curve(self.overonlylist,self.overonly_curve,
			self.param['overonly_lc_binsize'],pi_min=None,pi_max=None)

		self.overonly_curve_timeoffset = '%s/proc/%s_mgd_overonly_tim.flc' % (self.outdir,self.param['outbase'])
		cmd  = 'ftcalc infile=%s outfile=%s column=TIME expression="TIME+#TIMEZERO"\n' % (self.overonly_curve, self.overonly_curve_timeoffset)
		cmd += 'ftcalc infile=%s outfile=%s column=OVERONLY_RATE expression="RATE" clobber=yes\n' % (self.overonly_curve_timeoffset, self.overonly_curve_timeoffset)		
		cmd += 'ftcalc infile=%s outfile=%s column=OVERONLY_ERRIR expression="ERROR" clobber=yes\n' % (self.overonly_curve_timeoffset, self.overonly_curve_timeoffset)				
		cmd += 'ftcalc infile=%s outfile=%s column=OVERONLY_FRACEXP  expression="FRACEXP" clobber=yes\n' % (self.overonly_curve_timeoffset, self.overonly_curve_timeoffset)						
		cmd += 'fdelcol %s+1 colname="RATE" confirm=no proceed=yes\n'	% self.overonly_curve_timeoffset	
		cmd += 'fdelcol %s+1 colname="ERROR" confirm=no proceed=yes\n'	% self.overonly_curve_timeoffset	
		cmd += 'fdelcol %s+1 colname="FRACEXP" confirm=no proceed=yes\n'	% self.overonly_curve_timeoffset					
		print(cmd); os.system(cmd)
	"""

	"""
	def interporation_overonly2mkf(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'finterp %s+1 ' % self.merged_mkffile
		cmd += '%s+1 ' % self.overonly_curve_timeoffset
		cmd += '%s ' % self.merged_mkffile
		cmd += 'OVERONLY_RATE'
		print(cmd); os.system(cmd)			

		cmd  = 'fparkey %d ' % self.param['overonly_lc_binsize']
		cmd += '%s ' % self.merged_mkffile
		cmd += 'OVERBNSZ '
		cmd += 'comm="over-only light curve bin size (sec)" '
		cmd += 'add=yes '
		print(cmd); os.system(cmd)			
	"""

	def run_nimaketime(self,expr,outgti):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nimaketime '		
		cmd += 'infile=%s ' % self.merged_mkffile
		cmd += 'outfile=%s ' % outgti
		cmd += 'expr="%s" ' % expr 
		print(cmd); os.system(cmd)

	def run_nimaketime_nicersaa(self,outgti):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nimaketime '		
		cmd += 'infile=%s ' % self.merged_mkffile
		cmd += 'outfile=%s ' % outgti
		cmd += 'nicersaafilt=NO '
		cmd += 'expr="(NICER_SAA==1)" '
		print(cmd); os.system(cmd)		

	def merge_ufafiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.ufalist = '%s/proc/%s_ufalist.lis' % (self.outdir,self.param['outbase'])
		f = open(self.ufalist,'w')
		for niobs in self.niobs_list:
			sys.stdout.write('%s\n' % niobs.ufaevt_mpu7)
			f.write('%s\n' % niobs.ufaevt_mpu7)
		f.close()		

		self.merged_ufaevt = '%s/proc/%s_mgd_ufa.evt' % (self.outdir,self.param['outbase'])
		cmd = 'nimpumerge infiles=@%s outfile=%s mpulist=7' % (self.ufalist,self.merged_ufaevt)
		print(cmd); os.system(cmd)

	def run_nicerclean(self,ingti,outevt):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nicerclean infile=%s outfile=%s ' % (self.merged_ufaevt,outevt)
		if ingti != None:
			cmd += 'gtifile=%s ' % ingti
		print(cmd); os.system(cmd)

	"""
	def prepare_gtifiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.fgti_night = '%s/proc/gti/%s_night.gti' % (self.outdir,self.param['outbase'])
		self.run_nimaketime("(SUNSHINE==0)",self.fgti_night)

		self.fgti_day = '%s/proc/gti/%s_day.gti' % (self.outdir,self.param['outbase'])
		self.run_nimaketime("(SUNSHINE==1)",self.fgti_day)

		self.fgti_nicersaa = '%s/proc/gti/%s_nicersaa.gti' % (self.outdir,self.param['outbase'])
		self.run_nimaketime_nicersaa(self.fgti_nicersaa)		
		
		self.fgti_night_overcut_list = []
		self.fgti_day_overcut_list   = []		
		self.fgti_overcut_list   = []				
		for over_rates in self.param['overonly_rate_thresholds']:
			rate_min = over_rates[0]
			rate_max = over_rates[1]	

			index = self.param['overonly_rate_thresholds'].index(over_rates)		

			fgti_night = '%s/proc/gti/%s_night_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			self.fgti_night_overcut_list.append(fgti_night)
			expr = "(SUNSHINE==0)&&(OVERONLY_RATE>%.3f)&&(OVERONLY_RATE<=%.3f)" % (rate_min,rate_max)
			self.run_nimaketime(expr,fgti_night)

			fgti_day = '%s/proc/gti/%s_day_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			self.fgti_day_overcut_list.append(fgti_day)
			expr = "(SUNSHINE==1)&&(OVERONLY_RATE>%.3f)&&(OVERONLY_RATE<=%.3f)" % (rate_min,rate_max)
			self.run_nimaketime(expr,fgti_day)	

			fgti = '%s/proc/gti/%s_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			self.fgti_overcut_list.append(fgti)
			expr = "(OVERONLY_RATE>%.3f)&&(OVERONLY_RATE<=%.3f)" % (rate_min,rate_max)
			self.run_nimaketime(expr,fgti)						
	"""

	"""	
	def show_gtifiles_exposure(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		fgti_night_exp = get_total_gti_exposure(self.fgti_night,flag_dump=False,extension_name='STDGTI')
		fgti_day_exp = get_total_gti_exposure(self.fgti_day,flag_dump=False,extension_name='STDGTI')
		fgti_nicersaa_exp = get_total_gti_exposure(self.fgti_nicersaa,flag_dump=False,extension_name='STDGTI')		
		dump  = "----- Exposure ----- \n"
		dump += "                     : Night (sec)   Day (sec)\n"
		for over_rates in self.param['overonly_rate_thresholds']:
			index = self.param['overonly_rate_thresholds'].index(over_rates)		
			rate_min = over_rates[0]
			rate_max = over_rates[1]	
			exp_night = get_total_gti_exposure(self.fgti_night_overcut_list[index],flag_dump=False,extension_name='STDGTI')		
			exp_day = get_total_gti_exposure(self.fgti_day_overcut_list[index],flag_dump=False,extension_name='STDGTI')					
			dump += "Band-%d (%.1f-%.1f cps): %.1f      %.1f\n" % (index, rate_min, rate_max, exp_night, exp_day)
		dump += "Total                : %.1f      %.1f\n" % (fgti_night_exp,fgti_day_exp)
		dump += "NICER SAA: %.1f (s)\n" % fgti_nicersaa_exp

		print(dump)
		self.fgti_result = '%s/proc/gti/gti_exposure.txt' % self.outdir
		f = open(self.fgti_result,'w')
		f.write(dump)
		f.close()
	"""

	def set_response_files(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		if os.getenv('NICER_RESP_PATH') != '':
			self.rmffile = '%s/%s' % (os.getenv('NICER_RESP_PATH'),self.param['rmffile'])
			self.arffile = '%s/%s' % (os.getenv('NICER_RESP_PATH'),self.param['arffile'])
		else:
			self.rmffile = None
			self.arffile = None 

	"""
	def extract_cleaned_events_spectra(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		dict_keywords = {"MRGSTART":self.date_obs,"MRGSTOP":self.date_end}

		dict_keywords["FILTRCND"] = "standard"
		self.merged_clevt = self.merged_ufaevt.replace('_ufa.evt','_cl.evt')
		self.run_nicerclean(None,self.merged_clevt)
		self.clpha_all = xselect_extract_spectrum(self.merged_clevt,
			outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)		

		if get_total_gti_exposure(self.fgti_night,flag_dump=False,extension_name='STDGTI') > 0.0:
			dict_keywords["FILTRCND"] = "night-all"
			self.clevt_night = self.fgti_night.replace('.gti','.evt')
			self.run_nicerclean(self.fgti_night,self.clevt_night)
			self.clpha_night = xselect_extract_spectrum(self.clevt_night,
				outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)	
		if get_total_gti_exposure(self.fgti_day,flag_dump=False,extension_name='STDGTI') > 0.0:
			dict_keywords["FILTRCND"] = "day-all"
			self.clevt_day = self.fgti_day.replace('.gti','.evt')
			self.run_nicerclean(self.fgti_day,self.clevt_day)
			self.clpha_day = xselect_extract_spectrum(self.clevt_day,
				outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)	
		if get_total_gti_exposure(self.fgti_nicersaa,flag_dump=False,extension_name='STDGTI') > 0.0:
			dict_keywords["FILTRCND"] = "nicersaa"
			self.evt_nicersaa = self.fgti_nicersaa.replace('.gti','.evt')
			self.run_nicerclean(self.fgti_nicersaa,self.evt_nicersaa)
			self.pha_nicersaa = xselect_extract_spectrum(self.evt_nicersaa,
				outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)

		for over_rates in self.param['overonly_rate_thresholds']:
			rate_min = over_rates[0]
			rate_max = over_rates[1]	
			index = self.param['overonly_rate_thresholds'].index(over_rates)		

			fgti = '%s/proc/gti/%s_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			if get_total_gti_exposure(fgti,flag_dump=False,extension_name='STDGTI') > 0.0:
				dict_keywords["FILTRCND"] = "%.1e<over<=%.1e(all)" % (rate_min,rate_max)
				clevt = fgti.replace('.gti','.evt')
				self.run_nicerclean(fgti,clevt)
				clpha = xselect_extract_spectrum(clevt,
					outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)

			fgti_night = '%s/proc/gti/%s_night_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			if get_total_gti_exposure(fgti_night,flag_dump=False,extension_name='STDGTI') > 0.0:
				dict_keywords["FILTRCND"] = "%.1e<over<=%.1e(night)" % (rate_min,rate_max)
				clevt_night = fgti_night.replace('.gti','.evt')
				self.run_nicerclean(fgti_night,clevt_night)
				clpha_night = xselect_extract_spectrum(clevt_night,
					outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)

			fgti_day = '%s/proc/gti/%s_day_overcut%d.gti' % (self.outdir,self.param['outbase'],index)
			if get_total_gti_exposure(fgti_day,flag_dump=False,extension_name='STDGTI') > 0.0:
				dict_keywords["FILTRCND"] = "%.1e<over<=%.1e(day)" % (rate_min,rate_max)
				clevt_day = fgti_day.replace('.gti','.evt')
				self.run_nicerclean(fgti_day,clevt_day)
				clpha_day = xselect_extract_spectrum(clevt_day,
					outdir=None,rmffile=self.rmffile,arffile=self.arffile,dict_keywords=dict_keywords)
	"""

	"""
	def plot_spectrum(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		for phafile in glob.glob('%s/proc/gti/*.pha' % self.outdir):
			hdu = pyfits.open(phafile)
			exposure = float(hdu['SPECTRUM'].header['EXPOSURE'])
			otitle = "%s (%.3f ks)" % (phafile, exposure/1000.0)
			title  = ""
			ftitle = "%s (%s...%s)" % (hdu['SPECTRUM'].header['FILTRCND'],
				hdu['SPECTRUM'].header['MRGSTART'],hdu['SPECTRUM'].header['MRGSTOP'])
			plot_xspec_spectrum(phafile,otitle=otitle,title=title,ftitle=ftitle)

		for phafile in glob.glob('%s/proc/*.pha' % self.outdir):
			hdu = pyfits.open(phafile)
			exposure = float(hdu['SPECTRUM'].header['EXPOSURE'])
			otitle = "%s (%.3f ks)" % (phafile, exposure/1000.0)
			title  = ""
			ftitle = "%s (%s...%s)" % (hdu['SPECTRUM'].header['FILTRCND'],
				hdu['SPECTRUM'].header['MRGSTART'],hdu['SPECTRUM'].header['MRGSTOP'])
			plot_xspec_spectrum(phafile,otitle=otitle,title=title,ftitle=ftitle)
	"""

	def run_barycentric_correction(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.clevt_bary_lst = '%s/timing/%s_clevt_bary.lis' % (self.outdir,self.param['outbase'])
		f = open(self.clevt_bary_lst,'w')
		for niobs in self.niobs_list:
			niobs.run_barycentric_correction(ra=self.param['ra'],dec=self.param['dec'],ephem=self.param['ephem'])
			sys.stdout.write('%s\n' % niobs.clevt_bary)
			f.write('%s\n' % niobs.clevt_bary)
		f.close()		

	"""
	def extract_lowbackground_data(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.fgti_lowbgd = '%s/timing/%s_lowbgd.gti' % (self.outdir,self.param['outbase']) 
		expr = "OVERONLY_RATE<=%.3f" % self.param['pulse_search_overonly_threshold']
		self.run_nimaketime(expr,self.fgti_lowbgd)	

		self.flag_lowbgd_data = False 
		if get_total_gti_exposure(self.fgti_lowbgd,flag_dump=False,extension_name='STDGTI') > 0.0:
			self.flag_lowbgd_data = True
			self.clevt_lowbgd = self.fgti_lowbgd.replace('.gti','.evt')
			self.run_nicerclean(self.fgti_lowbgd,self.clevt_lowbgd)
			self.clpha_lowbgd = xselect_extract_spectrum(self.clevt_lowbgd,
				outdir=None,rmffile=self.rmffile,arffile=self.arffile)	
			self.clevt_lowbgd_bary = barycentric_correction(self.clevt_lowbgd,self.merged_orbfile,
				ra=self.param['ra'],dec=self.param['dec'],ephem=self.param['ephem'])
			self.clevt_lowbgd_bary_esel = nicer_fselect_energy(self.clevt_lowbgd_bary,
				emin_keV=self.param['pulse_search_emin_keV'],emax_keV=self.param['pulse_search_emax_keV'])

			self.clflc_lowbgd = '%s_%sto%skeV_%ds.flc' % (os.path.splitext(self.clevt_lowbgd)[0],
				str(self.param['lc_emin_keV']).replace('.','p'),
				str(self.param['lc_emax_keV']).replace('.','p'),
				self.param['lc_binsize'])
			xselect_extract_curve(self.clevt_lowbgd,self.clflc_lowbgd,
				self.param['lc_binsize'],
				pi_min=KEV_TO_PI*self.param['lc_emin_keV'],
				pi_max=KEV_TO_PI*self.param['lc_emax_keV'])
			plot_curve(self.clflc_lowbgd)
	"""
	
	"""
	def plot_curve(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
		self.merged_clflc_esel = '%s_%sto%skeV_%ds.flc' % (os.path.splitext(self.merged_clevt)[0],
			str(self.param['lc_emin_keV']).replace('.','p'),
			str(self.param['lc_emax_keV']).replace('.','p'),
			self.param['lc_binsize'])
		xselect_extract_curve(self.merged_clevt,self.merged_clflc_esel,
			self.param['lc_binsize'],
			pi_min=KEV_TO_PI*self.param['lc_emin_keV'],
			pi_max=KEV_TO_PI*self.param['lc_emax_keV'])
		plot_curve(self.merged_clflc_esel)
	"""

	def save_setup(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		if self.flag_multimode:
			cmd = 'cp %s %s ' % (self.filename_obsid_path,self.outdir)
			print(cmd);os.system(cmd)
		else:
			f = open('%s/input.lst' % self.outdir,'w')
			f.write('%s\n' % self.obsid_path)
			f.close()

		fparam_out = '%s/%s' % (self.outdir, os.path.basename(self.fparam))
		f = open(fparam_out,'w')
		f.write(yaml.dump(self.param))
		f.close()

	def run(self):
		sys.stdout.write('--run--\n')

		# Initialization 
		self.show_input_parameters()
		self.set_obsid_path_list()
		self.show_obsid_path_list()
		self.make_output_directory()
		self.load_parameterfile()
		self.show_input_parameters()
		self.set_nicer_observations()
		self.set_response_files()
		# Main process 
		self.run_nicerl2()
		self.run_niprefilter2()
		self.merge_mkffiles()
		self.merge_orbfiles()
		self.merge_ufafiles()
		#self.extract_overonly_event()         # Teru's method 		
		#self.extract_overonly_curve()         # Teru's method 		
		#self.interporation_overonly2mkf()     # Teru's method 		
		#self.prepare_gtifiles()               # Teru's method 		
		#self.show_gtifiles_exposure()         # Teru's method 		
		#self.extract_cleaned_events_spectra() # Teru's method 		
		#self.plot_spectrum()                  
		#self.plot_curve()                    
		# 
		# Timing 
		self.run_barycentric_correction() # for individual ObsID
		#self.extract_lowbackground_data() # for merged data 
		# End 
		self.save_setup()

	def run_org180319(self):
		sys.stdout.write('--run--\n')

		# Initialization 
		self.show_input_parameters()
		self.set_obsid_path_list()
		self.show_obsid_path_list()
		self.make_output_directory()
		self.load_parameterfile()
		self.show_input_parameters()
		self.set_nicer_observations()
		self.set_response_files()
		# Main process 
		#self.run_nicerl2()
		#self.run_niprefilter2()
		#self.merge_mkffiles()
		#self.merge_orbfiles()
		#self.merge_ufafiles()
		#self.extract_overonly_event()         # Teru's method 		
		#self.extract_overonly_curve()         # Teru's method 		
		#self.interporation_overonly2mkf()     # Teru's method 		
		#self.prepare_gtifiles()               # Teru's method 		
		#self.show_gtifiles_exposure()         # Teru's method 		
		#self.extract_cleaned_events_spectra() # Teru's method 		
		#self.plot_spectrum()                  
		#self.plot_curve()                    
		# 
		# Timing 
		#self.run_barycentric_correction() # for individual ObsID
		#self.extract_lowbackground_data() # for merged data 
		# End 
		self.save_setup()


