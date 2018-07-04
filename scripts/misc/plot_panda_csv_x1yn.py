#!/usr/bin/env python

import os 
import sys 
import yaml 
import argparse
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

import matplotlib as mpl
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '0' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

help_message = """
(usage) %s csvfile yaml
""" % sys.argv[0]
parser = argparse.ArgumentParser(description=help_message)
parser.add_argument('csvfile',metavar='csvfile',type=str,help='Input CSV file') 
parser.add_argument('yamlfile',metavar='yamlfile',type=str,help='Input yaml file') 
args = parser.parse_args()

if not os.path.exists(args.csvfile):
	sys.stderr.write('error: file %s does not exist.' % args.csvfile)
	quit()
df = pd.read_csv(args.csvfile)

if not os.path.exists(args.yamlfile):
	sys.stderr.write('error: file %s does not exist.' % args.yamlfile)
	quit()
param = yaml.load(open(args.yamlfile))

num_of_ycol = len(param['ycolumns'])

fig, axes = plt.subplots(num_of_ycol,1,sharex=True,sharey=False,
	figsize=(param['panel_size'][0],param['panel_size'][1]*num_of_ycol))
for i in range(num_of_ycol):
	x = df[param['xcolumn']].values - param['xoffset']
	y = (df[param['ycolumns'][i]].values)/float(param['ynorms'][i])
	if param['yerrors'][i][0] == "None":
		axes[i].plot(x,y,"o",markersize=8,markerfacecolor="r",markeredgecolor="k")	
	else:
		yerr_min = (df[param['yerrors'][i][0]].values)/float(param['ynorms'][i])
		yerr_max = (df[param['yerrors'][i][0]].values)/float(param['ynorms'][i])
		axes[i].errorbar(x,y,yerr=[yerr_min,yerr_max],
			fmt="o",markersize=8,markerfacecolor="r",markeredgecolor="k",
			color="k")		

	if i == num_of_ycol - 1:
		plt.xlabel(param['xlabel'])
	axes[i].set_xlim(param['xranges'])
	axes[i].set_ylabel(param['ylabels'][i])
	axes[i].set_autoscaley_on(False)
	axes[i].set_ylim(param['yranges'][i])
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("sample.pdf")
