#!/usr/bin/env python
# Usage: ./bgw-epsilon.py eps0mat.h5 epsmat.h5

import argparse
import asetk.format.bgw as bgw
import numpy as np
import os


# Define command line parser
parser = argparse.ArgumentParser(
    description='Reads dielectric matrix from BerkeleyGW output\
                 and plots macroscopic dielectric function.')
parser.add_argument('--version', action='version', version='%(prog)s 12.05.2016')
parser.add_argument(
    '--eps0',
    metavar='FILENAME', 
    default='eps0mat.h5',
    help='Filename of static dielectric matrix from BGW')
parser.add_argument(
    '--eps',
    metavar='FILENAME', 
    default=None,
    help='Filename of frequency-dependent dielectric matrix from BGW')
parser.add_argument(
    '--plot',
    action='store_true',
    default=False,
    help='Plot macroscopic dielectric function versus frequency.')

args = parser.parse_args()

eps_inv = bgw.DielectricMatrix.from_hdf5_db(eps0=args.eps0, eps=args.eps)

eps_macro = eps_inv.macroscopic_epsilon()
print("Macroscopic static dielectric constant: epsilon(0) = {:.6f}".format(eps_macro[0]))
#print(r"$\varepsilon(\omega=0) = 1/\varepsilon^{-1}_{GG'}(\vec{q}=0,\omega=0)$")


if args.plot:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    #...

