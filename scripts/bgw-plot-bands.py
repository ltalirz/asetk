#!/usr/bin/env python
# Usage: ./bgw-plot-bands.py eqp.dat

import argparse
import asetk.format.bgw as bgw
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pickle


# Define command line parser
parser = argparse.ArgumentParser(
    description='Reads energy levels from BerkeleyGW output\
                 and plots (1d) band structure to bands.png.')
parser.add_argument('--version', action='version', version='%(prog)s 29.01.2015')
parser.add_argument(
    'source',
    metavar='FILENAME', 
    help='Source file containing energy levels from BGW')
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
parser.add_argument(
    '--pickle',
    action='store_true',
    default=False,
    help='Save python spectrum object using pickle module.')

args = parser.parse_args()
filename = args.source

window = args.window * 2

print("Reading data from {f}".format(f=filename))

# we try to guess the format
if not args.format:
    ext = os.path.splitext(args.source)[-1]

    if args.source == "eqp.dat":
        args.format = 'eqp'
    elif ext == '.log':
        args.format = 'log'
    else:
        print("Error: Please specify source format using --format")

if args.format == 'eqp':
    spectrum   = bgw.Spectrum.from_eqp(filename, mode=args.mode)
elif args.format == 'log':
    spectrum   = bgw.Spectrum.from_log(filename, mode=args.mode)
else:
    print("Error: Unrecognized format specification {}".format(args.format))

print(spectrum)

data = None
for spin, dispersion in zip(spectrum.spins, spectrum.dispersions):
    for i in range(len(dispersion.energylevels)):
        E = dispersion.energylevels[i].energies
        fermi = dispersion.energylevels[i].fermi

        k = np.array([ dispersion.kvectors[i][2] for l_ in range(len(E))])
        # QE likes to place kpoints at -0.5 instead of +0.5
        #k = [ kp if kp >= 0 else kp+1 for kp in k]

        #k *= np.pi
        #E -= fermi
        fermi = 0

        #plt.plot(k,E, 'ko', markersize=2.0)
        plt.plot(k,E, 'ko')
        #plt.ylim(fermi - window/2, fermi + window/2)
        #plt.ylim(fermi - window , fermi)
        plt.xlabel(r'k [$\frac{2\pi}{a}$]')
        plt.ylabel('E [eV]')

        if data is not None:
            data = np.concatenate( (data, np.array(list(zip(k,E)))), axis=0 )
        else:
            data = np.array(list(zip(k,E)))

np.savetxt('bands.dat', data)

#plt.show()

pngname='bands.png'
print("Saving band structure plot to {f}".format(f=pngname))
plt.savefig(pngname, transparent=True, dpi=150)


if args.pickle:
    import pickle
    fname = 'spectrum.p'
    print("Saving pickled spectrum {}".format(fname))
    pickle.dump(spectrum, open(fname, "wb"), protocol=2)

