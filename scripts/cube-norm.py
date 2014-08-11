#!/usr/bin/env python
# usage: norm.py <filename>
# calculates norm of cube file
from asetk.format.cube import Cube
import asetk.atomistic.constants as constants
import numpy as np
import argparse


# Define command line parser
parser = argparse.ArgumentParser(
    description='Calculates norm of Gaussian cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 11.8.2014')
parser.add_argument(
    'cube',
    metavar='FILENAME',
    help='The Gaussian cube file')

args = parser.parse_args()


c = Cube.from_file(args.cube, read_data = True)
norm = np.sum(c.data) * c.dv * (constants.Angstrom / constants.a0)**3
print("{} has norm {}".format(args.cube, norm))
