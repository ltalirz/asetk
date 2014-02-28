#!/usr/bin/env python
# Produce density of states from CP2K MoLog file
# Usage: cp2k-dos.py cp2k.out
import argparse
import atk.format.cp2k as cp2k
import numpy as np


# Define command line parser
parser = argparse.ArgumentParser(
    description='Produce density of states from CP2K MoLog file')
parser.add_argument('--version', action='version', version='%(prog)s 28.02.2014')
parser.add_argument(
    'out',
    metavar='WILDCARD', 
    help='CP2K output containing energy levels')
parser.add_argument(
    '--plot',
    metavar='BOOL', 
    default=True,
    help='Whether to plot DOS.')
parser.add_argument(
    '--tofile',
    metavar='FILENAME', 
    default=None,
    help='Write DOS to file (if specified).')
parser.add_argument(
    '--window',
    metavar='ENERGY',
    default=3,
    help='plot range [-window,window] around Fermi.')
parser.add_argument(
    '--delta',
    metavar='ENERGY',
    default=0.001,
    help='the energy grid spacing')
parser.add_argument(
    '--sigma',
    metavar='ENERGY',
    default=0.075,
    help='The sigma of the Gaussian broadening. FWHM = sigma*sqrt(8*ln(2))')
parser.add_argument(
    '--nsigma',
    default=10,
    metavar='INT',
    help='Reduce computational cost by specifying up to how many sigma the Gaussians should be evaluated')

args = parser.parse_args()

spectrum = cp2k.Spectrum.from_mo(args.out)
print(spectrum)

E, DOS = spectrum.dos(sigma = args.sigma, nsigma = args.nsigma, deltaE = args.delta)



# Write dos to file
if args.tofile:
    header = "DOS from {}, sigma={}, nsigma={}".format(args.out, sigma,nsigma)
    np.savetxt(args.tofile, np.array(E,DOS).T, header=header)

# Plot DOS
if args.plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.plot(DOS, E)

    fermi = spectrum.fermi
    plt.plot([0, np.max(DOS)], [fermi, fermi], '--k')

    energies = spectrum.energies
    plt.plot(np.zeros(len(energies)), energies, 'ro')

    plt.xticks([])
    plt.xlim([0, np.max(DOS)])
    plt.ylim(fermi-args.window, fermi+args.window)
    plt.ylabel('E [eV]')
    plt.savefig(args.out+'-dos.png', dpi=200, transparent=True)


