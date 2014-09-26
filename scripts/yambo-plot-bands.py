#!/usr/bin/env python
# Usage: ./yambo-plot-bands.py yambo.out QP 

import argparse
import asetk.format.yambo as yambo
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os


# Define command line parser
parser = argparse.ArgumentParser(
    description='Reads energy levels from Yambo output\
                 and plots (1d) band structure to bands.png.')
parser.add_argument('--version', action='version', version='%(prog)s 12.09.2014')
parser.add_argument(
    'source',
    metavar='FILENAME', 
    help='Source file containing energy levels from Yambo')
parser.add_argument(
    '--format',
    metavar='STRING', 
    default=None,
    help='May be "netcdf_db" or "output" or "qp" (optional)')
parser.add_argument(
    '--mode',
    default='QP',
    metavar='STRING', 
    help='May be "DFT" or "QP" (optional)')
parser.add_argument(
    '--window',
    default=5,
    metavar='FLOAT', 
    type=float,
    help='Will plot [Ef-window, Ef+window]')

args = parser.parse_args()
filename = args.source

window = args.window * 2

print "Reading data from {f}".format(f=filename)

# we try to guess the format
if not args.format:
    ext = os.path.splitext(args.source)[-1]

    if args.source == "o.qp":
        args.format = 'qp'
    elif ext == '.QP':
        args.format = 'netcdf_db'
    elif ext == '.out':
        args.format = 'output'
    else:
        print("Error: Please specify source forat using --format")
 
if args.format == 'output' or args.format == 'qp':
    spectrum   = yambo.Spectrum.from_output(filename, mode=args.mode)
elif args.format == 'netcdf_db':
    spectrum   = yambo.Spectrum.from_netcdf_db(filename, mode=args.mode)
else:
    print("Error: Unrecognized format specification {}".format(args.format))

print spectrum

data = None
for spin, dispersion in zip(spectrum.spins, spectrum.dispersions):
    for i in range(len(dispersion.energylevels)):
        E = dispersion.energylevels[i].energies
        fermi = dispersion.energylevels[i].fermi

        k = np.array([ dispersion.kvectors[i][0] for l_ in range(len(E))])
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
