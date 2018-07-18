#!/usr/bin/env python

import sys 

if len(sys.argv) != 3 :
	sys.stderr.write('%s flux d_kpc\n' % sys.argv[0])
	exit()

flux  = float(sys.argv[1])
d_kpc = float(sys.argv[2])
luminosity = 1.2e+32 * (flux / 1e-12) * d_kpc**2 

dump  = "flux: %.3e (erg/s/cm2)\n" % flux
dump += "distance : %.3e (kpc)\n" % d_kpc
dump += "luminosity: %.3e (erg/s)" % luminosity
print(dump)