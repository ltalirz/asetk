#!/usr/bin/env python
import numpy as np
import argparse
from asetk.format.qe import QECube

# Define command line parser
parser = argparse.ArgumentParser(
    description='Converts intermediate cube-like format of pp.x to Gaussian cube format.')
parser.add_argument('--version', action='version', version='%(prog)s 13.10.2016')
parser.add_argument(
    'cubes',
    nargs='+',
    metavar='FILENAMES',
    help='QE cube file(s) to be converted. Conversion simply adds .cube extension.')

args = parser.parse_args()

for fname in args.cubes:
    print("Reading QE cube file {}".format(fname))
    c = QECube.from_file(fname, read_data=True)

    outname = "{}.cube".format(fname)
    print("Writting Gaussian cube file {}".format(outname))
    c.write_cube_file(outname)
    print("")

print("Job done")
