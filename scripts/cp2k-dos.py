#!/usr/bin/env python
# Produce density of states from CP2K MoLog file
# Usage: cp2k-dos.py cp2k.out
import argparse
import asetk.format.cp2k as cp2k
import numpy as np
import os


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
    default=None,
    help='The sigma of the Gaussian broadening. Equivalent to setting \
          FWHM = sigma*sqrt(8*ln(2)) (kept for backward compatibility).')
parser.add_argument(
    '--FWHM',
    metavar='ENERGY',
    default=0.1,
    help='Full-width half-maximum of broadening function.')
parser.add_argument(
    '--bmethod',
    default='Gaussian',
    metavar='STRING',
    help='Method used for broadening: "Gaussian" or "Lorentzian"')
parser.add_argument(
    '--bepsilon',
    default=1e-3,
    metavar='WEIGHT',
    help='Reduce computational cost by specifying quantiles that may be neglected.')

args = parser.parse_args()

lfname, lfext = os.path.splitext(args.out)
if lfext == '.MOLog':
        spectrum = cp2k.Spectrum.from_mo(args.out)
else:
    spectrum = cp2k.Spectrum.from_output(args.out)
print(spectrum)

if args.sigma:
    FWHM = np.sqrt(8.0 * np.log(2.0)) * args.sigma

for e,s in zip(spectrum.energylevels, spectrum.spins):
    E, DOS = e.dos(bmethod=args.bmethod, FWHM=args.FWHM, bepsilon=args.bepsilon, delta_e = args.delta)

    # Write dos to file
    if args.tofile:
        header = "DOS from {}, sigma={}, nsigma={}, spin="\
                .format(args.out, args.sigma, args.nsigma, s)
        np.savetxt("spin_{}_{}".format(s, args.tofile), 
                   np.array(zip(E,DOS)).T, header=header)

    # Plot DOS
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plt.plot(DOS, E, label='spin {}'.format(s))

        fermi = e.fermi
        plt.plot([0, np.max(DOS)], [fermi, fermi], '--k')

        energies = spectrum.energies
        plt.plot(np.zeros(len(energies)), energies, 'ro') 

        plt.xticks([])
        plt.xlim([0, np.max(DOS)])
        plt.ylim(fermi-args.window, fermi+args.window)
        plt.ylabel('E [eV]')

if args.plot:
    plt.legend()
    plt.savefig("spin_{}_dos.png".format(s, args.out),
                dpi=200, transparent=True)


