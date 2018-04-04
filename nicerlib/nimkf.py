import matplotlib
matplotlib.use('TkAgg')

import os 
import sys 
import yaml
import glob 
import pyfits

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

from pyheasoft import * 
from niconst import * 

#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['font.size'] = 18
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.xmargin'] = '0' #'.05'
plt.rcParams['axes.ymargin'] = '.1'

def merge_pdf(outpdf_list,outpdf_merged):
	cmd = 'convert '
	for outpdf in outpdf_list:
		cmd += '%s ' % outpdf
	cmd += '%s ' % outpdf_merged
	print(cmd);os.system(cmd)	
	print("open %s" % outpdf_merged)
	return outpdf_merged

class NicerMKF():
	def __init__(self,mkffile,outdir):
		self.mkffile = mkffile
		self.outdir  = outdir		
		sys.stdout.write('mkffile: %s\n' % self.mkffile)
		sys.stdout.write('outdir: %s\n' % self.outdir)		

		self.basename = os.path.splitext(os.path.basename(self.mkffile).replace('.gz',''))[0]
		self.outpdf_list = []

	def mkdir(self,recreate=True):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		if recreate:
			cmd = 'rm -rf %s;mkdir -p %s' % (self.outdir,self.outdir)
		print(cmd);os.system(cmd)

	def fplot_background_histogram(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		colname_list = [
			["FPM_DOUBLE_COUNT",200,0.0,200.0,True],
			["FPM_OVERONLY_COUNT",200,0.0,200.0,True],			
			["FPM_UNDERONLY_COUNT",200,0.0,200.0,True],						
			["FPM_FT_COUNT",200,0.0,200.0,True],									
			["FPM_NOISE25_COUNT",200,0.0,200.0,True],		
			["FPM_RATIO_REJ_COUNT",200,0.0,200.0,True],		
			["FPM_XRAY_PI_0000_0025",200,0.0,200.0,True],		
			["FPM_XRAY_PI_0025_0200",200,0.0,200.0,True],		
			["FPM_XRAY_PI_0200_0800",200,0.0,200.0,True],		
			["FPM_XRAY_PI_0800_1200",200,0.0,200.0,True],		
			["FPM_XRAY_PI_1200_1500",200,0.0,200.0,True],		
			["FPM_XRAY_PI_1500_1700",200,0.0,200.0,True],	
			["COR_SAX",20,0.0,20.0,False],		
			["SUN_ANGLE",180,0.0,180.0,False],		
			["ANG_DIST",180,0.0,180.0,False]
			]
		outpdf_list = []
		for i in colname_list:
			outpdf = '%s/%s_%s_hist.pdf' % (self.outdir,self.basename,i[0])
			cmd  = 'fplot_histogram.py %s %s ' % (self.mkffile,outpdf)
			cmd += '%s %d %.3f %.3f ' % (i[0],i[1],i[2],i[3])
			if i[4] == True:
				cmd += '--logx '
			cmd += '--logy '
			print(cmd);os.system(cmd)
			outpdf_list.append(outpdf)

		merge_pdf(outpdf_list,'%s/%s_hist_merge.pdf' % (self.outdir, self.basename))

	def fplot_background_scatter(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)
	
		colname_list = [
			["SUN_ANGLE","FPM_XRAY_PI_0000_0025",False,True],
			["FPM_UNDERONLY_COUNT","FPM_XRAY_PI_0000_0025",True,True],
			["FPM_NOISE25_COUNT","FPM_XRAY_PI_0000_0025",False,True],			
			["FPM_FT_COUNT","FPM_XRAY_PI_0000_0025",False,True],                        
			["BR_EARTH","FPM_XRAY_PI_0000_0025",False,True],                                    
			["FPM_OVERONLY_COUNT","FPM_XRAY_PI_0800_1200",True,True],				
			["FPM_OVERONLY_COUNT","FPM_XRAY_PI_1200_1500",True,True],							
			["FPM_OVERONLY_COUNT","FPM_XRAY_PI_1500_1700",True,True],
			["TIME","FPM_XRAY_PI_1200_1500",False,True],
			["SAT_LAT","FPM_XRAY_PI_1200_1500",False,True],
			["SAT_LON","FPM_XRAY_PI_1200_1500",False,True],            
			["SAT_ALT","FPM_XRAY_PI_1200_1500",False,True],                        
			]	
		outpdf_list = []
		for i in colname_list:
			outpdf = '%s/%s_vs_%s_scat.pdf' % (self.outdir,i[0],i[1])
			cmd  = 'fplot_scatter.py %s %s ' % (self.mkffile,outpdf)
			cmd += '%s %s ' % (i[0],i[1])	
			if i[2]:
				cmd += '--logx '
			if i[3]:
				cmd += '--logy '				
			print(cmd);os.system(cmd)
			outpdf_list.append(outpdf)
		
		merge_pdf(outpdf_list,'%s/%s_scatter_merge.pdf' % (self.outdir, self.basename))

	def plot_detid_map(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		detid_expmap = sum(self.hdu['PREFILTER'].data['FPM_ON'].astype(int))
		detid_mpu_double_count = sum(self.hdu['PREFILTER'].data['MPU_DOUBLE_COUNT'].astype(np.uint64))
		detid_mpu_ft_count = sum(self.hdu['PREFILTER'].data['MPU_FT_COUNT'].astype(np.uint64))
		detid_mpu_noise25_count = sum(self.hdu['PREFILTER'].data['MPU_NOISE25_COUNT'].astype(np.uint64))
		detid_mpu_overonly_count = sum(self.hdu['PREFILTER'].data['MPU_OVERONLY_COUNT'].astype(np.uint64))
		detid_mpu_underonly_count = sum(self.hdu['PREFILTER'].data['MPU_UNDERONLY_COUNT'].astype(np.uint64))
		#detid_mpu_ft_pi_ave = sum(self.hdu['PREFILTER'].data['MPU_FT_PI_AVG'].astype(np.float128))
		#detid_mpu_ft_pi_err = sum(self.hdu['PREFILTER'].data['MPU_FT_PI_ERR'].astype(np.float128))
		#detid_mpu_ft_pi_fast_ave = sum(self.hdu['PREFILTER'].data['MPU_FT_PI_FAST_AVG'].astype(np.float128))
		#detid_mpu_ft_pi_fast_err = sum(self.hdu['PREFILTER'].data['MPU_FT_PI_FAST_ERR'].astype(np.float128))
		#detid_mpu_noise25_pi_avg = sum(self.hdu['PREFILTER'].data['MPU_NOISE25_PI_AVG'].astype(np.float128))
		#detid_mpu_noise25_pi_err = sum(self.hdu['PREFILTER'].data['MPU_NOISE25_PI_ERR'].astype(np.float128))

		detid_mpu_double_rate = detid_mpu_double_count.astype(np.float128) / detid_expmap.astype(np.float128)
		detid_mpu_ft_rate = detid_mpu_ft_count.astype(np.float128) / detid_expmap.astype(np.float128)
		detid_mpu_noise25_rate = detid_mpu_noise25_count.astype(np.float128) / detid_expmap.astype(np.float128)
		detid_mpu_overonly_rate = detid_mpu_overonly_count.astype(np.float128) / detid_expmap.astype(np.float128)
		detid_mpu_underonly_rate = detid_mpu_underonly_count.astype(np.float128) / detid_expmap.astype(np.float128)
		#detid_mpu_ft_pi_ave_timeave = detid_mpu_ft_pi_ave.astype(np.float128) / detid_expmap.astype(np.float128)
		#detid_mpu_ft_pi_err_timeave = detid_mpu_ft_pi_err.astype(np.float128) / detid_expmap.astype(np.float128)		
		#detid_mpu_ft_pi_fast_ave_timeave = detid_mpu_ft_pi_fast_ave.astype(np.float128) / detid_expmap.astype(np.float128)		
		#detid_mpu_ft_pi_fast_err_timeave = detid_mpu_ft_pi_fast_err.astype(np.float128) / detid_expmap.astype(np.float128)		
		#detid_mpu_noise25_pi_avg_timeave = detid_mpu_noise25_pi_avg.astype(np.float128) / detid_expmap.astype(np.float128)		
		#detid_mpu_noise25_pi_err_timeave = detid_mpu_noise25_pi_err.astype(np.float128) / detid_expmap.astype(np.float128)										

		detid_map_list = [
			[detid_expmap,'Exposure Map'],
			[detid_mpu_double_rate,'MPU DOUBLE Average Rate'],
			[detid_mpu_ft_rate,'MPU FT Average Rate'],
			[detid_mpu_noise25_rate,'MPU Noise25 Average Rate'],
			[detid_mpu_overonly_rate,'MPU OVERONLY Average Rate'],
			[detid_mpu_underonly_rate,'MPU UNDERONLY Average Rate']
			#[detid_mpu_ft_pi_ave,'MPU FT PI AVG time-average'],	
			#[detid_mpu_ft_pi_err_timeave,'MPU FT PI ERR time-average'],
			#[detid_mpu_ft_pi_fast_ave_timeave,'MPU FT PI FAST AVG time-average'],
			#[detid_mpu_ft_pi_fast_err_timeave,'MPU FT PI FAST ERR time-average'],
			#[detid_mpu_noise25_pi_avg_timeave,'MPU NOISE25 PI AVG time-average'],
			#[detid_mpu_noise25_pi_err_timeave,'MPU NOISE25 PI AVG time-average']
			]
		for detid_map in detid_map_list:
			fig = plt.figure()
			sns_plot = sns.heatmap(detid_map[0],
				annot=True,annot_kws={"size":8,"color":'r'},fmt='g',cmap='Blues')
			sns_plot.set_xlabel("FPM (DET_ID=10xMPU+1xFPM)")
			sns_plot.set_ylabel("MPU")
			sns_plot.set_title(detid_map[1])
			pdfname = '%s/%s_detmap.pdf' % (self.outdir,detid_map[1].replace(' ',''))
			fig.savefig(pdfname)

			self.outpdf_list.append(pdfname)

	def set_panda_dataframe(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.hdu = pyfits.open(self.mkffile)

		self.dframe = pd.DataFrame({
			'SAT_LAT':self.hdu['PREFILTER'].data['SAT_LAT'].astype('float'),	
			'SAT_LON':self.hdu['PREFILTER'].data['SAT_LON'].astype('float'),	
			'SAT_ALT':self.hdu['PREFILTER'].data['SAT_ALT'].astype('float'),		
			'BR_EARTH':self.hdu['PREFILTER'].data['BR_EARTH'].astype('float'),	
			'SUNSHINE':self.hdu['PREFILTER'].data['SUNSHINE'].astype('int'),	
			'SUNSHINE_STR':self.hdu['PREFILTER'].data['SUNSHINE'].astype('string'),				
			'SUN_ANGLE':self.hdu['PREFILTER'].data['SUN_ANGLE'].astype('float'),	
			'MOON_ANGLE':self.hdu['PREFILTER'].data['MOON_ANGLE'].astype('float'),	
			'ANG_DIST':self.hdu['PREFILTER'].data['ANG_DIST'].astype('float'),	
			'SAA_TIME':self.hdu['PREFILTER'].data['SAA_TIME'].astype('float'),	
			'COR_SAX':self.hdu['PREFILTER'].data['COR_SAX'].astype('float'),						
			'NICER_SAA':self.hdu['PREFILTER'].data['NICER_SAA'].astype('int'),
			'NICER_SAA_STR':self.hdu['PREFILTER'].data['NICER_SAA'].astype('string'),			
			'FPM_RATIO_REJ_COUNT':self.hdu['PREFILTER'].data['FPM_RATIO_REJ_COUNT'].astype('float'),	
			'FPM_XRAY_PI_0000_0025':self.hdu['PREFILTER'].data['FPM_XRAY_PI_0000_0025'].astype('float'),
			'FPM_XRAY_PI_0025_0200':self.hdu['PREFILTER'].data['FPM_XRAY_PI_0025_0200'].astype('float'),
			'FPM_XRAY_PI_0200_0800':self.hdu['PREFILTER'].data['FPM_XRAY_PI_0200_0800'].astype('float'),
			'FPM_XRAY_PI_0800_1200':self.hdu['PREFILTER'].data['FPM_XRAY_PI_0800_1200'].astype('float'),
			'FPM_XRAY_PI_1200_1500':self.hdu['PREFILTER'].data['FPM_XRAY_PI_1200_1500'].astype('float'),
			'FPM_XRAY_PI_1500_1700':self.hdu['PREFILTER'].data['FPM_XRAY_PI_1500_1700'].astype('float'),					
			'FPM_DOUBLE_COUNT':self.hdu['PREFILTER'].data['FPM_DOUBLE_COUNT'].astype('float'),						
			'FPM_OVERONLY_COUNT':self.hdu['PREFILTER'].data['FPM_OVERONLY_COUNT'].astype('float'),						
			'FPM_UNDERONLY_COUNT':self.hdu['PREFILTER'].data['FPM_UNDERONLY_COUNT'].astype('float'),						
			'FPM_FT_COUNT':self.hdu['PREFILTER'].data['FPM_FT_COUNT'].astype('float'),						
			'FPM_NOISE25_COUNT':self.hdu['PREFILTER'].data['FPM_NOISE25_COUNT'].astype('float')
			})
		print(self.dframe)

	def plot_parameter_correlation_optical_loading(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		sns.set(rc={'figure.figsize':(11.7,8.27)})

		sns_plot_optical = sns.PairGrid(self.dframe, 
			vars=["FPM_UNDERONLY_COUNT","FPM_NOISE25_COUNT","BR_EARTH","SUN_ANGLE","MOON_ANGLE"],
			hue="SUNSHINE_STR")
		sns_plot_optical.map(plt.scatter)
		plt.title("%s optical loading " % os.path.basename(self.mkffile))				
		pngname = '%s/%s_optical_scat.png' % (self.outdir,self.basename)
		#pdfname = '%s/%s_optical_scat.pdf' % (self.outdir,self.basename)		
		sns_plot_optical.savefig(pngname)
		#sns_plot_optical.savefig(pdfname)

		#self.outpdf_list.append(pdfname)
		self.outpdf_list.append(pngname)

	def plot_parameter_correlation_particle_background(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		sns_plot_particle = sns.PairGrid(self.dframe, 
			vars=["FPM_OVERONLY_COUNT","FPM_XRAY_PI_0800_1200","FPM_XRAY_PI_1200_1500","FPM_RATIO_REJ_COUNT","COR_SAX","SAT_LON"],
			hue="NICER_SAA_STR")
		sns_plot_particle.map(plt.scatter)
		plt.title("%s scattering " % os.path.basename(self.mkffile))		
		pngname = '%s/%s_praticlebkg_scat.png' % (self.outdir,self.basename)
		#pdfname = '%s/%s_praticlebkg_scat.pdf' % (self.outdir,self.basename)		
		sns_plot_particle.savefig(pngname)	
		#sns_plot_particle.savefig(pdfname)			

		self.outpdf_list.append(pngname)		

	def plot_jointplot(self,xcol,ycol,outpdfname,title=""):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		sns.set(rc={'figure.figsize':(11.7,8.27)})

		outpdf_path = '%s/%s' % (self.outdir,outpdfname)
		sns_plot = sns.jointplot(x=xcol,y=ycol,data=self.dframe)
		plt.title(title)
		sns_plot.savefig(outpdf_path)

		self.outpdf_list.append(outpdf_path)		

	def plot_correlations(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.plot_jointplot(xcol="FPM_OVERONLY_COUNT",ycol="FPM_XRAY_PI_0800_1200",
			outpdfname="FPM_OVERONLY_COUNT_vs_FPM_XRAY_PI_0800_1200.pdf",
			title="%s" % os.path.basename(self.mkffile))

		self.plot_jointplot(xcol="FPM_UNDERONLY_COUNT",ycol="FPM_NOISE25_COUNT",
			outpdfname="FPM_UNDERONLY_COUNT_vs_FPM_NOISE25_COUNT.pdf",
			title="%s" % os.path.basename(self.mkffile))

	def merge_pdffiles(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		self.outpdf_merge = '%s/%s_mkf_merge.pdf' % (self.outdir,self.basename)
		cmd  = 'convert '
		for outpdf in self.outpdf_list:
			cmd += '%s ' % outpdf 
		cmd += '%s ' % self.outpdf_merge
		print(cmd);os.system(cmd)

	def test(self):
		sys.stdout.write('=== %s ===\n' % sys._getframe().f_code.co_name)

		sns.set(rc={'figure.figsize':(11.7,8.27)})

		sns_plot_optical = sns.PairGrid(self.dframe, 
			vars=["FPM_UNDERONLY_COUNT","FPM_NOISE25_COUNT","BR_EARTH","SUN_ANGLE","MOON_ANGLE"],
			hue="SUNSHINE_STR")
		sns_plot_optical.map(plt.scatter)
		plt.title("%s optical loading " % os.path.basename(self.mkffile))				
		pngname = '%s/%s_optical_scat.png' % (self.outdir,self.basename)
		sns_plot_optical.savefig(pngname)

	def run(self):
		self.mkdir()
		self.set_panda_dataframe()
		#self.fplot_background_histogram()
		#self.fplot_background_scatter()
		self.plot_detid_map()
		self.plot_parameter_correlation_optical_loading()
		self.plot_parameter_correlation_particle_background()
		self.plot_correlations()
		self.merge_pdffiles()
		#self.test()

"""
	  Column Name                Format     Dims       Units     TLMIN  TLMAX
	  1 TIME                       1D                  s
	  2 POSITION                   3E                  km
	  3 VELOCITY                   3E                  km/s
	  4 QUATERNION                 4E
	  5 PNTUNIT                    3E
	  6 POLAR                      3E                  rad, rad, km
	  7 RA                         1E                  deg
	  8 DEC                        1E                  deg
	  9 ROLL                       1E                  deg
	 10 SAT_LAT                    1E                  deg
	 11 SAT_LON                    1E                  deg
	 12 SAT_ALT                    1E                  km
	 13 ELV                        1E                  deg
	 14 BR_EARTH                   1E                  deg
	 15 SUNSHINE                   1I
	 16 FOV_FLAG                   1I
	 17 SUN_ANGLE                  1E                  deg
	 18 MOON_ANGLE                 1E                  deg
	 19 RAM_ANGLE                  1E                  deg
	 20 ANG_DIST                   1E                  deg
	 21 SAA                        1I
	 22 SAA_TIME                   1E                  s
	 23 COR_ASCA                   1E                  GeV/c
	 24 COR_SAX                    1E                  GeV/c
	 25 MCILWAIN_L                 1E
	 26 SUN_RA                     1E                  deg
	 27 SUN_DEC                    1E                  deg
	 28 MOON_RA                    1E                  deg
	 29 MOON_DEC                   1E                  deg
	 30 EARTH_RA                   1E                  deg
	 31 EARTH_DEC                  1E                  deg
	 32 TIME_ADJ                   1D                  s
	 33 ST_BBO                     1B
	 34 ST_VALID                   1B
	 35 ST_OBJECTS                 1B
	 36 ST_VIDEO_VDC               1B                  V
	 37 ATT_ANG_AZ                 1D                  deg
	 38 ATT_ANG_EL                 1D                  deg
	 39 RA_CMD                     1D                  deg
	 40 DEC_CMD                    1D                  deg
	 41 ATT_ERR_AZ                 1D                  deg
	 42 ATT_ERR_EL                 1D                  deg
	 43 ATT_STATE                  1B
	 44 ATT_MODE                   1B
	 45 ATT_SUBMODE_AZ             1B
	 46 ATT_SUBMODE_EL             1B
	 47 TARG_CMD                   1I
	 48 PPS_SOURCE                 1J
	 49 PPS_ERR_LOWPASS            1D                  s
	 50 GPS_INIT                   1B
	 51 GPS_CONVERGED              1B
	 52 NICER_SAA                  1B
	 53 ST_STARS                   1B
	 54 ST_FAILCODE                1B
	 55 MPU_ALL_COUNT              56I     (8,7)
	 56 MPU_OVER_COUNT             56I     (8,7)
	 57 MPU_UNDER_COUNT            56I     (8,7)
	 58 MPU_XRAY_COUNT             56I     (8,7)
	 59 TOT_ALL_COUNT              1J
	 60 TOT_UNDER_COUNT            1J
	 61 TOT_OVER_COUNT             1J
	 62 TOT_XRAY_COUNT             1J
	 63 FPM_ON                     56L     (8,7)
	 64 NUM_FPM_ON                 1J
	 65 FPM_RATIO_REJ_COUNT        1E
	 66 FPM_XRAY_PI_0000_0025      1E
	 67 FPM_XRAY_PI_0025_0200      1E
	 68 FPM_XRAY_PI_0200_0800      1E
	 69 FPM_XRAY_PI_0800_1200      1E
	 70 FPM_XRAY_PI_1200_1500      1E
	 71 FPM_XRAY_PI_1500_1700      1E
	 72 MPU_DOUBLE_COUNT           56J     (8,7)
	 73 MPU_FT_COUNT               56I     (8,7)
	 74 MPU_NOISE25_COUNT          56I     (8,7)
	 75 MPU_OVERONLY_COUNT         56J     (8,7)
	 76 MPU_UNDERONLY_COUNT        56J     (8,7)
	 77 FPM_DOUBLE_COUNT           1E
	 78 FPM_OVERONLY_COUNT         1E
	 79 FPM_UNDERONLY_COUNT        1E
	 80 FPM_FT_COUNT               1E
	 81 FPM_NOISE25_COUNT          1E
	 82 MPU_FT_PI_AVG              56E     (8,7)       chan
	 83 MPU_FT_PI_ERR              56E     (8,7)       chan
	 84 MPU_FT_PI_FAST_AVG         56E     (8,7)       chan
	 85 MPU_FT_PI_FAST_ERR         56E     (8,7)       chan
	 86 MPU_NOISE25_PI_AVG         56E     (8,7)       chan
	 87 MPU_NOISE25_PI_ERR         56E     (8,7)       chan
"""












