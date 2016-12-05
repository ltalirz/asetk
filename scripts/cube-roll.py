#!/usr/bin/env python
import numpy as np
import argparse
from asetk.format.cube import Cube

# Define command line parser
parser = argparse.ArgumentParser(
    description='Rolls cube file by specified shift along direction.')
parser.add_argument('--version', action='version', version='%(prog)s 12.05.2016')
parser.add_argument(
    'cube',
    metavar='filename',
    help='Cube file to be sliced.')
parser.add_argument(
    '--dir',
    metavar='direction',
    type=str,
    help='"x", "y" or "z"')
parser.add_argument(
    '--shift',
    metavar='INTEGER',
    type=int,
    default=None,
    help='Shift specified in steps of cube file grid.')
parser.add_argument(
    '--distance',
    metavar='FLOAT',
    type=float,
    help='Shift specified in Angstroms.')

args = parser.parse_args()

# Make list of jobs
print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

c.roll(dir=args.dir, shift=args.shift, distance=args.distance)

fname_out = 's.' + args.cube
print("Writing {}".format(fname_out))
c.write_cube_file(fname_out)
