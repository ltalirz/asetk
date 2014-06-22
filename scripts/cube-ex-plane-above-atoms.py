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
    '--height',
    metavar='DISTANCE',
    default=2.5,
    type=float,
    help='The height [Angstroms] above the topmost atom, where to extract the plane.')

args = parser.parse_args()

# Make list of jobs
print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane = c.get_plane_above_atoms(args.height)

outfile = '{f}.z{i}.dat'.format(f=args.cube, i=args.height)
print("Writing plane data to {}".format(outfile))
np.savetxt(outfile, plane)
