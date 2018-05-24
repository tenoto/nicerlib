import os 
import sys 
import yaml
import glob 
#import pyfits
import numpy as np 
import astropy.io.fits as pyfits

from pyheasoft import * 
from niconst import * 
from nibgd import * 

NICER_DATA_SUBDIRECTORY_LIST = ['auxil','log','xti'] 
NICER_OUTPUT_SUBDIRS = ['proc','timing']
NICER_YAML_KEYWORDS = ["title","outbase","ra","dec","expr_str","expr","rmffile","arffile","flag_nicerl2","flag_niprefilter2"]

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
			if os.path.exists(self.catfile.replace('.gz','')):
				cmd += 'gzip -f %s\n' % (self.catfile.replace('.gz',''))						
			print(cmd);os.system(cmd)

			self.clevt_mpu7  = '%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (self.obsid_path, self.obsid)
			self.ufaevt_mpu7 = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt.gz' % (self.obsid_path, self.obsid)		

			cmd = 'rm -f %s/xti/event_cl/ni%s_0mpu[0-6]_ufa.evt' % (self.obsid_path,self.obsid)
			print(cmd);os.system(cmd)

	def run_niprefilter2(self):
		sys.stdout.write('=== %s (ObsID=%s) ===\n' % (sys._getframe().f_code.co_name,self.obsid))

		cmd  = 'niprefilter2 '
		cmd += 'indir=%s ' % self.obsid_path
		cmd += 'infile=%s ' % self.mkffile
		cmd += 'outfile=%s  ' % self.mkffile.replace('.gz','')
		cmd += 'clobber=yes '
		print(cmd);os.system(cmd)

		cmd  = 'rm -f %s; gzip %s\n' % (self.mkffile,self.mkffile.replace('.gz',''))
		print(cmd);os.system(cmd)

	def make_output_directory(self,outdir,flag_recreate=True):
		sys.stdout.write('=== %s (ObsID=%s) ===\n' % (sys._getframe().f_code.co_name,self.obsid))

		self.outdir = outdir
		if flag_recreate:
			cmd = 'rm -rf %s;mkdir -p %s' % (self.outdir,self.outdir)
		else:
			cmd = 'mkdir -p %s' % (self.outdir,self.outdir)			
		print(cmd);os.system(cmd)

	def nimaketime(self,expr,outgti,options=""):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nimaketime '		
		cmd += 'infile=%s ' % self.mkffile
		cmd += 'outfile=%s ' % outgti
		cmd += 'expr="%s" ' % expr 
		cmd += options 
		print(cmd); os.system(cmd)

		self.outgti = outgti 
		return outgti 

	def nicerclean(self,ingti,outevt):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nicerclean infile=%s outfile=%s ' % (self.ufaevt_mpu7,outevt)
		if ingti != None:
			cmd += 'gtifile=%s ' % ingti
		print(cmd); os.system(cmd)

		self.outevt = outevt
		return outevt

	def run_barycentric_correction(self,ra,dec,ephem=""):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.clevt_bary = barycentric_correction(self.clevt_mpu7,self.orbfile,ra=ra,dec=dec,ephem=ephem)

	def photonphase(self,inevt,outbaryevt,tempo_parfile):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'cp %s %s' % (inevt,outbaryevt)
		print(cmd);os.system(cmd)

		#cmd  = 'photonphase --addphase --barytime '
		cmd  = 'photonphase --addphase '
		#cmd += '--ephem %s ' % self.param['ephem']
		cmd += '--orbfile %s ' % self.orbfile
		cmd += '--plotfile %s ' % outbaryevt.replace('.evt','.png')
		cmd += '%s %s >& ' % (outbaryevt,tempo_parfile)
		cmd += '%s ' % outbaryevt.replace('.evt','.log')
		print(cmd);os.system(cmd)

		"""
		self.photonphase_evt = '%s/ni%s_%s_baryphase.evt' % (self.outdir,self.obsid,self.param['expr_str'])
		self.photonphase_log = '%s/ni%s_%s_baryphase.log' % (self.outdir,self.obsid,self.param['expr_str'])
		self.photonphase_png = '%s/ni%s_%s_baryphase.png' % (self.outdir,self.obsid,self.param['expr_str'])		

		cmd  = 'cp %s %s' % (self.cl2evt_mpu7,self.photonphase_evt)
		print(cmd);os.system(cmd)

		cmd  = 'photonphase --addphase '
		cmd += '--ephem %s ' % self.param['ephem']
		cmd += '--orbfile %s ' % self.orbfile
		cmd += '--plotfile %s ' % self.photonphase_png
		cmd += '%s %s >& ' % (self.photonphase_evt,self.param['parfile'])
		cmd += '%s ' % self.photonphase_log
		print(cmd);os.system(cmd)
		"""

	def make_background_spectrum(self,inevt,outdir):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.nibgd = NicerBackgroundModel(self.obsid_path,outdir)
		self.nibgd.set_cleaned_event(inevt)
		self.nibgd.run()

		exit()

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
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)	
		for key in NICER_YAML_KEYWORDS:
			if not self.param.has_key(key):
				sys.stdout.write('yaml file does not have a keyword %s.\n' % key)
				quit()

	def set_nicer_observations(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)		
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

	def merge_ufafiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		#self.ufalist = '%s/proc/%s_ufalist.lis' % (self.outdir,self.param['outbase'])
		self.merged_ufalst = '%s/proc/%s_mgd_ufa.lst' % (self.outdir,self.param['outbase'])
		f = open(self.merged_ufalst,'w')
		for niobs in self.niobs_list:
			sys.stdout.write('%s\n' % niobs.ufaevt_mpu7)
			f.write('%s\n' % niobs.ufaevt_mpu7)
		f.close()		

		self.merged_ufaevt = '%s/proc/%s_mgd_ufa.evt' % (self.outdir,self.param['outbase'])
		cmd = 'nimpumerge infiles=@%s outfile=%s mpulist=7' % (self.merged_ufalst,self.merged_ufaevt)
		print(cmd); os.system(cmd)

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

	def run_nicerclean(self,ingti,outevt):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		cmd  = 'nicerclean infile=%s outfile=%s ' % (self.merged_ufaevt,outevt)
		if ingti != None:
			cmd += 'gtifile=%s ' % ingti
		print(cmd); os.system(cmd)

	def run_each_obsid(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		subdir = '%s/obsid' % self.outdir
		cmd  = 'rm -rf %s; mkdir -p %s;' % (subdir,subdir)
		print(cmd); os.system(cmd)

		for niobs in self.niobs_list:
			targetdir = '%s/%s' % (subdir,niobs.obsid)
			niobs.make_output_directory(targetdir,flag_recreate=True)
			outgti = '%s/ni%s_%s.gti' % (targetdir,niobs.obsid,self.param['expr_str'])
			niobs.nimaketime(expr=self.param['expr'],outgti=outgti,options="")
			outevt = '%s/ni%s_0mpu7_%s.evt' % (targetdir,niobs.obsid,self.param['expr_str'])
			niobs.nicerclean(ingti=outgti,outevt=outevt)
			if self.param['flag_photonphase']:
				outbaryevt = outevt.replace('.evt','_bary.evt')
				if not os.path.exists(self.param['tempo_parfile']):
					sys.stderr.write('file %s does not exist.\n' % self.param['tempo_parfile'])
					quit()
				niobs.photonphase(outevt,outbaryevt,tempo_parfile=self.param['tempo_parfile'])
			if self.param['flag_bgdspectrum']:
				bgddir = '%s/bgd' % targetdir
				niobs.make_background_spectrum(inevt=outevt,outdir=bgddir)

	def set_response_files(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		if os.getenv('NICER_RESP_PATH') != '':
			self.rmffile = '%s/%s' % (os.getenv('NICER_RESP_PATH'),self.param['rmffile'])
			self.arffile = '%s/%s' % (os.getenv('NICER_RESP_PATH'),self.param['arffile'])
		else:
			self.rmffile = None
			self.arffile = None 

	def extract_spectrum(self,evtfile):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		outpha = evtfile.replace('.evt','.pha')
		cmd  = 'fxselect_extract_spec.py '
		cmd += '--inputevtfits=%s ' % evtfile
		cmd += '--outputpha=%s ' % outpha
		cmd += '--rmf=%s ' % self.rmffile
		cmd += '--arf=%s ' % self.arffile
		print(cmd);os.system(cmd)

		return outpha


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

		self.show_input_parameters()
		self.set_obsid_path_list()
		self.show_obsid_path_list()
		self.make_output_directory()
		self.load_parameterfile()
		self.show_input_parameters()
		self.set_nicer_observations()
		self.set_response_files()
		# Main process 
		if self.param['flag_nicerl2']:
			self.run_nicerl2()
		if self.param['flag_niprefilter2']:
			self.run_niprefilter2()
		if self.param['flag_merge_first']:
			self.merge_mkffiles()
			self.merge_orbfiles()
			self.merge_ufafiles()		
			self.gtifile = '%s/proc/%s_mgd_%s.gti' % (self.outdir,self.outbase,self.param['expr_str'])
			self.merged_cl2evt = self.gtifile.replace('.gti','.evt').replace('_ufa','_cl2')
			self.run_nimaketime(expr=self.param['expr'],outgti=self.gtifile)
			self.run_nicerclean(ingti=self.gtifile,outevt=self.merged_cl2evt)
			self.merged_cl2pha = self.extract_spectrum(self.merged_cl2evt)
		else:
			self.run_each_obsid()
		self.save_setup()


