#!/usr/bin/env python
# Reads energy levels from qe output and plots (1d) band structure
# to file "bands.png"

import asetk.format.qe as qe
from sys import argv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
import argparse

# Define command line parser
parser = argparse.ArgumentParser(
    description='Extract and plot (1d) band structure from QE')
parser.add_argument('--version', action='version', version='%(prog)s 24.03.2014')
parser.add_argument(
    'prefix',
    metavar='STRING', 
    help='Prefix of .save directory to be read.')
parser.add_argument(
    '--plot',
    action='store_true',
    default=True,
    help='Whether to plot bands to bands.png.')
parser.add_argument(
    '--window',
    metavar='ENERGY',
    default=3,
    type=float,
    help='Plot range [-window,window] around Fermi.')
parser.add_argument(
    '--save_bands',
    action='store_true',
    default=False,
    help='Save each band to file b<index>.dat')
parser.add_argument(
    '--pickle',
    action='store_true',
    default=False,
    help='Save python spectrum object using pickle module.')

args = parser.parse_args()

prefix = args.prefix

spectrum = qe.Spectrum.from_save(args.prefix)
dispersion = spectrum.dispersions[0]

data = None
k    = None
klist = []


for i in range(dispersion.nkpt):
    kpt = dispersion.kpoints[i]
    lev = kpt.energylevels

    E = lev.energies
    fermi = lev.fermi
    print(fermi)


    #k = np.array([ dispersion.kvectors[i][0] for l_ in range(len(E))])
    #k = [ kp if kp >= 0 else kp+1 for kp in k]
    # QE likes to place kpoints at -0.5 instead of +0.5
    if kpt.kvector[0] < 0:
        kpt.kvector[0] += 1.0

    if k is None:
        k = [0 for l_ in range(len(E))]
        klist.append(0)
    else:
        d = np.linalg.norm(dispersion.kvectors[i] - dispersion.kvectors[i-1])
        k += np.array([ d for l_ in range(len(E)) ] )
        klist.append(klist[i-1] + d)

    #k *= np.pi
    E -= fermi
    fermi = 0

    plt.plot(k,E, 'ko')
    window = [fermi -args.window, fermi+ args.window]
    plt.ylim(window)
    #plt.xlabel(r'k [$\frac{1}{\AA}$]')
    plt.xlabel(r'k [$\frac{2\pi}{a}$]')
    plt.ylabel('E [eV]')

    datablock = [ np.concatenate( (kpt.kvector, [e]) ) for e in E]
    if data is not None:
        data = np.concatenate( (data, datablock), axis=0 )
    else:
        data = datablock

np.savetxt(prefix+'_bands.dat', data, header='#kx  ky  kz  E[eV]', fmt='%.4e %.4e %.4e %.6e')
    
# alternative: grouping by band
if args.save_bands:
    bands = np.array([kpt.energies for kpt in dispersion.kpoints])
    for i in range(len(bands.T)):
        line = bands.T[i]
        fname = "{}_b{:03d}.dat".format(prefix, i+1)
        np.savetxt(fname,  zip(klist, line), header="k   E [eV]")

if args.pickle:
    import pickle
    fname = prefix + '_spectrum.p'
    pickle.dump(spectrum, open(fname, "wb"), protocol=2)

#plt.show()


plt.savefig(prefix+'_bands.png', transparent=True, dpi=150)

