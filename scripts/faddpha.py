#!/usr/bin/env python
# add pha files 

import sys
import os 
import math

"""
cd input 
rm -f ni10201301_o1u20cut_clscr_gtisel_merge.pha
mathpha << EOF
@src_pha.lst
C
ni10201301_o1u20cut_clscr_gtisel_merge.pha
CALC
NULL
0
EOF
"""

__author__  = "Teru Enoto"
__version__ = "1.0.0"
__email__   = "teruaki.enoto@gmail.com"
__status__  = "test"

if len(sys.argv) != 2:
	sys.stderr.write("usage: %s input_pha.lst\n" % sys.argv[0])
	quit()

fname_inlst = sys.argv[1]
inpha_lst = []
for line in open(fname_inlst):
	inpha_lst.append(line.split()[0])

fname_expr = 'tmp_expr.txt'
cmd = 'rm -f %s' % fname_expr
print(cmd);os.system(cmd)

f = open(fname_expr,'w')
for i in range(len(inpha_lst)):
	cmd = 'ln -s %s tmp_%s' % (inpha_lst[i],os.path.basename(inpha_lst[i]))
	print(cmd);os.system(cmd)
	if i < len(inpha_lst)-1:
		f.write("'tmp_%s' + " % os.path.basename(inpha_lst[i]))
	else:
		f.write("'tmp_%s'" % os.path.basename(inpha_lst[i]))			
f.close()	

fname_outpha = '%s_merge.pha' % os.path.splitext(sys.argv[1])[0]

cmd  = 'mathpha <<EOF\n'
cmd += '@%s\n' % fname_expr
cmd += 'C\n'
cmd += '%s\n' % fname_outpha
cmd += 'CALC\n'
cmd += 'NULL\n'
cmd += '0\n'
cmd += 'EOF\n'
print(cmd);os.system(cmd)

cmd = 'rm -f %s' % fname_expr
print(cmd);os.system(cmd)
for i in range(len(inpha_lst)):
	cmd = 'rm -f tmp_%s' % (os.path.basename(inpha_lst[i]))
	print(cmd);os.system(cmd)


"""
inqdp = sys.argv[2]
if not os.path.exists(inpha):
	sys.stderr.write("input pha file %s does not exists.\n" % inpha)
	quit()
if not os.path.exists(inqdp):
	sys.stderr.write("input qdp file %s does not exists.\n" % inqdp)
	quit()
if os.path.splitext(inqdp)[-1] not in ['.qdp']:
	sys.stderr.write("input qdp file must be qdp file (with .qdp extension): %s\n" % inqdp) 
	quit()
outpha = '%s.pha' % os.path.splitext(inqdp)[0]
if os.path.exists(outpha):
	sys.stderr.write("output pha file has already existed: %s\n" % outpha)
	quit()

cmd  = 'grppha<<EOF\n'
cmd += '%s\n' % inpha
cmd += '%s\n' % outpha

flagBody = False
for line in open(inqdp):
	cols = line.split()
	if cols[0] == '!':
		flagBody = True
		continue
	if not flagBody:
		continue
	start   = int(float(cols[0]) - float(cols[1]))
	stop    = int(float(cols[0]) + float(cols[1]) - 1)
	binsize = int(2.0*float(cols[1]))
	cmd += 'group %d %d %d\n' % (start, stop, binsize)
	if flagSystematic:
		binsyserr = math.sqrt(float(binsize)) * syserr
		cmd += 'systematics %d-%d %.4f \n' % (start,stop,binsyserr)

cmd += "exit\n"	
cmd += "EOF\n"	
os.system(cmd)

f = open('temp_header.txt','w')
f.write('HISTORY ---------------\n')
if flagSystematic:
	f.write('HISTORY fgrppha.py %s %s %s\n' % (sys.argv[1],sys.argv[2],sys.argv[3]))
else:	
	f.write('HISTORY fgrppha.py %s %s\n' % (sys.argv[1],sys.argv[2]))
for line in cmd.split('\n'):
	dump = 'HISTORY %s\n' % line
	f.write(dump)
f.close()

cmd = ''
for i in range(2):
	cmd += 'fthedit %s+%d \@temp_header.txt\n' % (outpha,i)
	if flagSystematic:
		cmd += 'fparkey %.4f "%s[%d]" SYSERR comm="systematic error for the spectral bins." add=yes\n' % (
			syserr, outpha, i)
cmd += 'rm -f temp_header.txt'	
os.system(cmd)	
"""
