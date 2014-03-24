#!/usr/bin/env python
# Usage: ./qe-exbands.py qe.out
# Reads energy levels from qe output and plots (1d) band structure
# to file "bands.png"
# Version 24.03.2014

import atk.format.qe as qe
from sys import argv
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt

prefix = argv[1]
window = 4 # window around fermi in eV

spectrum = qe.Spectrum.from_save(prefix)
dispersion = spectrum.dispersions[0]

data = None
k    = None
for i in range(len(dispersion.kpoints)):
    kpt = dispersion.kpoints[i]
    E = kpt.energies
    fermi = kpt.fermi


    #k = np.array([ dispersion.kvectors[i][0] for l_ in range(len(E))])
    # QE likes to place kpoints at -0.5 instead of +0.5
    #k = [ kp if kp >= 0 else kp+1 for kp in k]
    if k is None:
        k = [0 for l_ in range(len(E))]
    else:
        d = np.linalg.norm(dispersion.kvectors[i] - dispersion.kvectors[i-1])
        k += np.array([ d for l_ in range(len(E)) ] )

    #k *= np.pi
    E -= fermi
    fermi = 0

    plt.plot(k,E, 'ko')
    plt.ylim(fermi - window/2, fermi + window/2)
    #plt.xlabel(r'k [$\frac{1}{\AA}$]')
    plt.xlabel(r'k [$\frac{2\pi}{a}$]')
    plt.ylabel('E [eV]')

    datablock = [ np.concatenate( (dispersion.kvectors[i], [e]) ) for e in E]
    if data is not None:
        data = np.concatenate( (data, datablock), axis=0 )
    else:
        data = datablock

np.savetxt('bands.dat', data, header='#kx  ky  kz  E[eV]', fmt='%.4e %.4e %.4e %.6e')
    
#plt.show()

plt.savefig('bands.png', transparent=True, dpi=150)

