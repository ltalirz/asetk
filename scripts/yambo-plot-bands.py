#!/usr/bin/env python
# Usage: ./yambo-plot-bands.py yambo.out QP 

import argparse
import asetk.format.yambo as yambo
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


# Define command line parser
parser = argparse.ArgumentParser(
    description='Reads energy levels from Yambo output\
                 and plots (1d) band structure to bands.png.')
parser.add_argument('--version', action='version', version='%(prog)s 27.01.2014')
parser.add_argument(
    'out',
    metavar='FILENAME', 
    help='Yambo output containing energy levels')
parser.add_argument(
    'mode',
    metavar='STRING', 
    help='May be "DFT" or "QP" (optional)')

args = parser.parse_args()
filename = args.out

window = 10 # window around fermi in eV

print "Reading data from {f}".format(f=filename)
spectrum   = yambo.Spectrum.from_output(filename, mode=args.mode)

print spectrum

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
