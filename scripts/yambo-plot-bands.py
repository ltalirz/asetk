#!/usr/bin/env python
# Usage: ./yambo-plot-bands.py yambo.out
# Reads energy levels from Yambo output and plots (1d) band structure
# to file "bands.png"
# Version 27.1.2014

import atk.format.yambo as yambo
from sys import argv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

filename = argv[1]
window = 10 # window around fermi in eV

print "Reading data from {f}".format(f=filename)
spectrum   = yambo.Spectrum.from_output(filename)

data = None
for spin, dispersion in zip(spectrum.spins, spectrum.dispersions):
    for i in range(len(dispersion.energylevels)):
        E = dispersion.energylevels[i].energies
        fermi = dispersion.energylevels[i].fermi

        k = np.array([ dispersion.kvectors[i][2] for l_ in range(len(E))])
        # QE likes to place kpoints at -0.5 instead of +0.5
        k = [ kp if kp >= 0 else kp+1 for kp in k]

        #k *= np.pi
        #E -= fermi
        fermi = 0

        #plt.plot(k,E, 'ko', markersize=2.0)
        plt.plot(k,E, 'ko')
        plt.ylim(fermi - window/2, fermi + window/2)
        #plt.ylim(fermi - window , fermi)
        plt.xlabel(r'k [$\frac{2\pi}{a}$]')
        plt.ylabel('E [eV]')

        if data is not None:
            data = np.concatenate( (data, np.array(zip(k,E))), axis=0 )
        else:
            data = np.array(zip(k,E))

np.savetxt('bands.dat', data)

#plt.show()

pngname='bands.png'
print "Saving band structure plot to {f}".format(f=pngname)
plt.savefig(pngname, transparent=True, dpi=150)
