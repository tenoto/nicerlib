#!/usr/bin/env python

import os 

cmd  = 'rm -f ni1013010110_clmerge_screen_gtisel_bgdsub.pha'
print(cmd);os.system(cmd)

cmd  = "mathpha properr=yes errmeth='POISS-3' <<EOF\n" 
cmd += "ni1013010110_clmerge_screen_gtisel.pha-ni1013010110_mitbgd_BGMod_3C50.pha\n"
cmd += "R\n"
cmd += "ni1013010110_clmerge_screen_gtisel_bgdsub.pha\n"
cmd += "ni1013010110_clmerge_screen_gtisel.pha\n"
cmd += "%\n"
cmd += "0\n"
cmd += "EOF\n"
print(cmd);os.system(cmd)