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