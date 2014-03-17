#!/usr/bin/env python
import numpy as np
import argparse
from atk.format.cube import Cube

# Define command line parser
parser = argparse.ArgumentParser(
    description='Extracts plane from cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 27.01.2014')
parser.add_argument(
    'cube',
    metavar='',
    help='Cube file to be sliced.')
parser.add_argument(
    'dir',
    metavar='direction',
    type=str,
    help='Plane normal. May be "x", "y" or "z"')
parser.add_argument(
    'index',
    metavar='',
    type=int,
    help='Plane index.')

args = parser.parse_args()

# Make list of jobs
print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane = c.get_plane(dir=args.dir, i=args.index)

outfile = '{f}.plane{i}.dat'.format(f=args.cube, i=args.index)
print("Writing plane data to {}".format(outfile))
np.savetxt(outfile, plane)
