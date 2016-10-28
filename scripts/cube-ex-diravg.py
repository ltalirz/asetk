#!/usr/bin/env python
import numpy as np
import argparse
from asetk.format.cube import Cube

# Define command line parser
parser = argparse.ArgumentParser(
    description='Extracts directional average of cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 27.01.2014')
parser.add_argument(
    'cube',
    metavar='filename',
    help='Cube file to be sliced.')
parser.add_argument(
    'dir',
    metavar='direction',
    type=str,
    help='Plane normal. May be "x", "y" or "z"')

args = parser.parse_args()

# Make list of jobs
print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane = c.get_avg(dir=args.dir)

outfile = '{f}.diravg{d}.dat'.format(f=args.cube, d=args.dir)
print("Writing plane data to {}".format(outfile))
np.savetxt(outfile, plane)
