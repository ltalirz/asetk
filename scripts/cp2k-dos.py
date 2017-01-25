#!/usr/bin/env python
# Produce density of states from CP2K MOLog file
# Usage: cp2k-dos.py cp2k.out
import argparse
import asetk.format.cp2k as cp2k
import numpy as np
import os


# Define command line parser
parser = argparse.ArgumentParser(
    description='Produce density of states from CP2K output or .MOLog file')
parser.add_argument('--version', action='version', version='%(prog)s 24.01.2017')
parser.add_argument(
    'out',
    metavar='WILDCARD', 
    help='CP2K output containing energy levels (or: .MOLog / .pdos file)')
parser.add_argument(
    '--plot',
    metavar='BOOL', 
    default=True,
    help='Whether to plot DOS.')
parser.add_argument(
    '--to_file',
    dest='to_file',
    action='store_true',
    default=False,
    help='Whether to write DOS to file.')
parser.add_argument(
    '--window',
    metavar='ENERGY',
    default=3,
    type=float,
    help='plot range [-window,window] around Fermi [eV].')
parser.add_argument(
    '--delta',
    metavar='ENERGY',
    default=0.001,
    type=float,
    help='the energy grid spacing [eV]')
parser.add_argument(
    '--sigma',
    metavar='ENERGY',
    default=None,
    type=float,
    help='The sigma of the Gaussian broadening. Equivalent to setting \
          FWHM = sigma*sqrt(8*ln(2)) (kept for backward compatibility) [eV].')
parser.add_argument(
    '--FWHM',
    metavar='ENERGY',
    type=float,
    default=0.1,
    help='Full-width half-maximum of broadening function [eV].')
parser.add_argument(
    '--bmethod',
    default='Gaussian',
    metavar='STRING',
    help='Method used for broadening: "Gaussian" or "Lorentzian"')
parser.add_argument(
    '--bepsilon',
    type=float,
    default=1e-3,
    metavar='WEIGHT',
    help='Reduce computational cost by specifying quantiles that may be neglected.')

args = parser.parse_args()

lfname, lfext = os.path.splitext(args.out)
if lfext == '.MOLog':
        spectrum = cp2k.Spectrum.from_mo(args.out)
elif lfext == '.pdos':
    spectrum = cp2k.Spectrum.from_pdos(args.out)
else:
    spectrum = cp2k.Spectrum.from_output(args.out)
print(spectrum)

if args.sigma:
    args.FWHM = np.sqrt(8.0 * np.log(2.0)) * args.sigma

for e,s in zip(spectrum.energylevels, spectrum.spins):
    E, DOS = e.dos(bmethod=args.bmethod, FWHM=args.FWHM, bepsilon=args.bepsilon, delta_e = args.delta)

    dosint = np.sum(DOS) * args.delta
    print("Integrated DOS spin {}: {:.3f} electrons".format(s+1,dosint))
    dosmax = np.max(DOS)
    print("Maximum DOS spin {}: {:.3e} electrons / eV".format(s+1,dosmax))

    # Write dos to file
    if args.to_file:
        header = "DOS from {}, FWHM={} eV, spin={}\n E [eV]             DOS"\
                .format(args.out, args.FWHM, s+1)
        np.savetxt("dos_spin{}.dat".format(s+1), 
                   np.array(zip(E,DOS)), header=header)

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
        plt.xlim([0, dosmax])
        plt.ylim(fermi-args.window, fermi+args.window)
        plt.ylabel('E [eV]')

if args.plot:
    plt.legend()
    plt.savefig("dos.png", dpi=300, transparent=True)


