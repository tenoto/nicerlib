#!/usr/bin/env python

import os 
import sys 
import yaml 
import argparse
import subprocess
import numpy as np 
import pandas as pd 
import astropy.io.fits as pyfits 

COLUMNS_INPUT = ['data_id','src_pha','bgd_pha',
	'rmffile','arffile','model_xcm',
	'bin_minsig','bin_maxbin','fit_emin','fit_emax','param_free']
COLUMNS_ADD = ['exposure','total_rate',
	'tstart_met','tstop_met',
	'tstart_mjd_utc','tstop_mjd_utc','mjd_utc','mjd_width',
	'tstart_iso8601_utc','tstop_iso8601_utc',
	'grp_pha','rate_fitband','rate_error_fitband',
	'rchi2','chi2','dof','prob']	

def get_timeconv(met,outtype='MJD_UTC',outcol=2):
	cmd = 'rm -f tmp.log\n'
	print(cmd);os.system(cmd)

	cmd = 'nitimeconv.py %.6f -f met -s met > tmp.log' % met
	print(cmd);os.system(cmd)

	for line in open('tmp.log'):
		cols = line.split()
		if len(cols) == 0:
			continue 
		if cols[0] == outtype:
			outime = cols[outcol]			

	cmd = 'rm -f tmp.log\n'
	print(cmd);os.system(cmd)

	return outime

def get_MJD_UTC(met):
	return float(get_timeconv(met,outtype='MJD_UTC',outcol=2))

def get_ISO8601_UTC(met):
	return get_timeconv(met,outtype='ISO8601_UTC:',outcol=1)

def get_fitrange_rate(logfile):
	for line in open(logfile):
		cols = line.split()
		if "#Net count rate" in line:
			rate = float(cols[6])
			rate_err = float(cols[8])
			break
	return rate, rate_err

def get_chisqure(logfile):
	for line in open(logfile):
		cols = line.split()
		if "#Test statistic : Chi-Squared =" in line:
			chi2 = float(cols[5])
		if "# Reduced chi-squared =" in line:
			print(line)
			rchi2 = float(cols[4])
			dof = int(cols[6])
		if "# Null hypothesis probability =" in line:
			print(line)
			prob = float(cols[5])
	return rchi2, chi2, dof, prob

def get_flux(logfile,emin,emax):
	flag_found = False
	for line in open(logfile):
		cols = line.split() 
		if flag_found:
			flux_min = float(cols[6].replace('(',''))
			flux_max = float(cols[8].replace(')',''))
			break 
		if not "Model Flux" in line:
			continue
		emin_data = float(cols[8].replace('(',''))
		emax_data = float(cols[10])
		if emin == emin_data and emax == emax_data:
			print(line)
			flux = float(cols[5].replace('(',''))
			flag_found = True 
	flux_err_min = flux_min - flux
	flux_err_max = flux_max - flux 
	if not flag_found:
		print("no corresponding flux...")
		flux = np.nan
		flux_err_min = np.nan
		flux_err_max = np.nan 
	return flux, flux_err_min, flux_err_max

def get_parameter(logfile,npar):
	flag_show_pa = False
	for line in open(logfile):
		cols = line.split()
		if flag_show_pa and len(cols) > 1:
			if cols[1] == str(npar):
				value = float(cols[-3])
				break
		if "#Parameters defined:" in line:
			flag_show_pa = True
	if not flag_show_pa:
		value = np.nan

	flag_err = False
	for line in open(logfile):
		cols = line.split()
		if flag_err and len(cols) == 5:
			if cols[1] == str(npar) and cols[4][0] == '(' and cols[4][-1] == ')':
				print(line)
				tmp1,tmp2 = cols[4].split(',')
				err_min = float(tmp1.replace('(',''))
				err_max = float(tmp2.replace(')',''))
				break
		if "# Parameter   Confidence Range (1)" in line:
			flag_err = True
	if not flag_err:
		err_min = np.nan
		err_max = np.nan		
	return value, err_min, err_max

class XspecDataSet():
	def __init__(self,dataset,param):
		print("========== generate a new XspecDataSet ==========")
		self.dataset = dataset
		self.param = param 
		print(self.dataset)

		self.result = []
		for colname in COLUMNS_INPUT:
			self.result.append(self.dataset[colname])

	def bin_spec(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		basename = '%s_bin%ds%d' % (
			os.path.splitext(os.path.basename(self.dataset['src_pha']))[0],
			self.dataset['bin_minsig'],
			self.dataset['bin_maxbin'])
		grp_qdp = '%s.qdp' % basename
		self.grp_pha = '%s.pha' % basename

		cmd  = 'xspec <<EOF\n'
		cmd += 'data 1 %s\n' % self.dataset['src_pha']
		if not pd.isnull(self.dataset['bgd_pha']):
			cmd += 'back 1 %s\n' % self.dataset['bgd_pha']
		cmd += 'setplot energy\n'
		cmd += 'ignore 1:**-%.2f %.2f-**\n' % (self.dataset['fit_emin'],self.dataset['fit_emax'])
		cmd += 'setplot rebin %d %d\n' % (self.dataset['bin_minsig'],self.dataset['bin_maxbin'])
		cmd += 'setplot channel\n'
		cmd += 'iplot data\n'
		cmd += ''
		cmd += 'we %s\n' % basename
		cmd += 'exit\n'
		cmd += 'exit\n'
		print(cmd);os.system(cmd)

		cmd  = 'fgrppha.py %s %s' % (self.dataset['src_pha'],grp_qdp)
		print(cmd);os.system(cmd)

		cmd  = 'rm -f %s/%s' % (os.path.dirname(self.dataset['src_pha']),self.grp_pha)
		print(cmd);os.system(cmd)

		cmd  = 'mv %s %s;' % (self.grp_pha,os.path.dirname(self.dataset['src_pha']))
		cmd += 'rm -f %s.{qdp,pco};' % basename
		print(cmd);os.system(cmd)

		self.grp_pha = '%s/%s' % (os.path.dirname(self.dataset['src_pha']),self.grp_pha)
		self.result.append(self.grp_pha)

	def make_pcofile(self):
		self.fpco = self.grp_pha.replace('.pha','_fit.pco')
		f = open(self.fpco,'w')
		dump  = 'time off\n'
		dump += 'lwid 5 \n'
		dump += 'lwid 5 on 1..100 \n'	
		dump += 'csize 1.1\n'
		dump += 'lab pos y 2.8\n'
		dump += 'r x %.1f %.1f\n' % (float(self.param['plot_xmin']),float(self.param['plot_xmax']))
		dump += 'r y %.1e %.1e\n' % (float(self.param['plot_ymin']),float(self.param['plot_ymax']))
		dump += 'r y2 %.1f %.1f\n' % (float(self.param['plot_y2min']),float(self.param['plot_y2max']))
		dump += 'col 2 on 2\n'
		dump += 'win 2\n'
		dump += 'LAB  2 COL 2 LIN 0 100 JUS Lef POS 0.200000003 0 " "\n'
		f.write(dump)
		f.close()

	def check_header_keywords(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		hdu = pyfits.open(self.dataset['src_pha'])

		self.total_count = np.sum(hdu['SPECTRUM'].data['COUNTS'])
		self.exposure   = float(hdu['SPECTRUM'].header['EXPOSURE'])
		self.total_rate = float(self.total_count)/float(self.exposure)
		self.tstart_met = float(hdu['GTI'].data[0]['START'])
		self.tstop_met = float(hdu['GTI'].data[-1]['STOP'])

		self.tstart_mjd_utc = get_MJD_UTC(self.tstart_met)
		self.tstop_mjd_utc = get_MJD_UTC(self.tstop_met)		

		self.mjd_utc = 0.5*(self.tstart_mjd_utc+self.tstop_mjd_utc)
		self.mjd_width = 0.5*(-self.tstart_mjd_utc+self.tstop_mjd_utc)

		self.tstart_iso8601_utc = get_ISO8601_UTC(self.tstart_met)
		self.tstop_iso8601_utc = get_ISO8601_UTC(self.tstop_met)		

		self.result.append(self.exposure)
		self.result.append(self.total_rate)
		self.result.append(self.tstart_met)
		self.result.append(self.tstop_met)		
		self.result.append(self.tstart_mjd_utc)
		self.result.append(self.tstop_mjd_utc)
		self.result.append(self.mjd_utc)		
		self.result.append(self.mjd_width)
		self.result.append(self.tstart_iso8601_utc)				
		self.result.append(self.tstop_iso8601_utc)	

		self.obsid = hdu['SPECTRUM'].header['OBS_ID']
		self.objectname = hdu['SPECTRUM'].header['OBJECT']

		self.title = '%s %s %s (%.1f s)' % (
			self.objectname, self.obsid, self.tstart_iso8601_utc, self.exposure)
		self.subtitle = '%.3f cps (total), %.1f-%.1f keV fitted by %s' % (
			self.total_rate, 
			self.dataset['fit_emin'],self.dataset['fit_emax'],
			os.path.basename(self.dataset['model_xcm']))

	def fit_spec(self):
		sys.stdout.write('----- %s -----\n' % sys._getframe().f_code.co_name)

		self.flog_fit = self.grp_pha.replace('.pha','_fit.log')
		fxcm_fit = self.grp_pha.replace('.pha','_fit.xcm')		
		fps_fit = self.grp_pha.replace('.pha','_fit.ps')		

		cmd = 'rm -f %s %s %s' % (self.flog_fit,fxcm_fit,fps_fit)
		print(cmd);os.system(cmd)

		number_of_param_error = []
		for i in self.dataset['param_free'].split(','):
			number_of_param_error.append(int(i))

		cmd  = 'xspec<<EOF\n'
		cmd += 'data 1 %s\n' % self.grp_pha
		if not pd.isnull(self.dataset['bgd_pha']):
			cmd += 'back 1 %s\n' % self.dataset['bgd_pha']
		cmd += 'setplot energy\n'
		cmd += 'ignore 1:**-%0.3f %0.3f-**\n' % (self.dataset['fit_emin'],self.dataset['fit_emax'])
		cmd += '@%s\n' % self.dataset['model_xcm']
		cmd += 'query yes\n'	
		cmd += 'fit\n'		
		cmd += 'log %s\n' % self.flog_fit 
		cmd += 'save all %s\n' % fxcm_fit				
		cmd += 'show rate\n'
		cmd += 'show pa\n'		
		cmd += 'show fit\n'
		for eband in self.param['fit_flux_eband_list']:
			cmd += 'flux %.1f %.1f err 100 68.3\n' % (eband[0],eband[1])
		for n_err in number_of_param_error:
			cmd += 'err 1.0 %d\n' % n_err			
		cmd += 'log none\n'
		cmd += 'iplot ld del\n'
		cmd += '@%s\n' % self.fpco
		cmd += 'win 1 \n'
		cmd += 'la t %s\n' % self.title
		cmd += 'la f %s\n' % self.subtitle
		cmd += 'hard %s/cps\n' % fps_fit 
		cmd += 'exit\n'
		cmd += 'exit\n'
		cmd += 'EOF\n'
		print(cmd);os.system(cmd)

		cmd = 'ps2pdf %s' % fps_fit
		print(cmd);os.system(cmd)
		fpdf_fit = os.path.basename(fps_fit).replace('.ps','.pdf')
		cmd = 'mv %s %s; rm -f %s' % (fpdf_fit,os.path.dirname(fps_fit),fps_fit)
		print(cmd);os.system(cmd)

	def get_fit_parameters(self):
		self.rate_fitband, self.rate_error_fitband = get_fitrange_rate(self.flog_fit)
		self.result.append(self.rate_fitband)
		self.result.append(self.rate_error_fitband)

		self.rchi2, self.chi2, self.dof, self.prob = get_chisqure(self.flog_fit)
		self.result.append(self.rchi2)
		self.result.append(self.chi2)
		self.result.append(self.dof)
		self.result.append(self.prob)

		for eband in self.param['fit_flux_eband_list']:
			flux, flux_err_min, flux_err_max = get_flux(self.flog_fit,eband[0],eband[1])
			self.result.append(flux)
			self.result.append(flux_err_min)
			self.result.append(flux_err_max)
		for num in self.param['param_num']:
			value, err_min, err_max = get_parameter(self.flog_fit,num)
			self.result.append(value)
			self.result.append(err_min)	
			self.result.append(err_max)	

	def run(self):
		self.check_header_keywords()		
		self.bin_spec()
		self.make_pcofile()
		self.fit_spec()
		self.get_fit_parameters()

class FittingManager():
	def __init__(self,csvfile,yamlfile):
		self.csvfile = csvfile 
		self.yamlfile = yamlfile 

		if not os.path.exists(self.csvfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.csvfile)
			exit()			
		if not os.path.exists(self.yamlfile):
			sys.stderr.write('error: file %s does not exist.\n' % self.yamlfile)
			exit()						
		print("csvfile: %s" % self.csvfile)
		print("yamlfile: %s" % self.yamlfile)

		self.df = pd.read_csv(self.csvfile)
		self.param = yaml.load(open(self.yamlfile))

		self.outcsvfile = os.path.basename(self.csvfile).replace('.csv','_fit.csv')
		print("out csvfile: %s" % self.outcsvfile)

		self.COLUMNS_ADD2 = []
		for eband in self.param['fit_flux_eband_list']:
			colname = 'flux%sto%s' % (
				str(eband[0]).replace('.','p'),
				str(eband[1]).replace('.','p'))
			self.COLUMNS_ADD2.append(colname)
			self.COLUMNS_ADD2.append(colname+'_err_min')
			self.COLUMNS_ADD2.append(colname+'_err_max')			

		if len(self.param['param_num']) != len(self.param['param_name']):
			sys.stderr.write('number mismatch param_num and param_name\n')
			exit()
		for param_name in self.param['param_name']:
			self.COLUMNS_ADD2.append(param_name)
			self.COLUMNS_ADD2.append(param_name+'_err_min')			
			self.COLUMNS_ADD2.append(param_name+'_err_max')						
		print(self.COLUMNS_ADD2)

	def run(self):		

		self.csv_data_list = []
		for index, dataset in self.df.iterrows():
			dset = XspecDataSet(dataset,self.param)
			dset.run()
			self.csv_data_list.append(dset.result)
		print(self.csv_data_list)

		COLUMNS_ALL = COLUMNS_INPUT + COLUMNS_ADD + self.COLUMNS_ADD2

		df_out = pd.DataFrame(self.csv_data_list,columns=COLUMNS_ALL)
		df_out.to_csv(self.outcsvfile)			

		print(len(COLUMNS_ALL))
		print(len(self.csv_data_list[0]))

if __name__ == '__main__':
	cmd = 'which fhelp'
	resp = subprocess.Popen(cmd.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	if len(resp.stdout.readlines()) == 0:
		sys.stderr.write('error: set HEASOFT (e.g., heainit)\n')
		exit()

	help_message = """
	(example) %s csvfile
	""" % sys.argv[0]
	parser = argparse.ArgumentParser(description=help_message)
	parser.add_argument('csvfile',metavar='csvfile',type=str,help='Input pha csvfile.') 
	parser.add_argument('yamlfile',metavar='yamlfile',type=str,help='Input parameter yamlfile.') 
	args = parser.parse_args()

	fit = FittingManager(args.csvfile,args.yamlfile)
	fit.run()
