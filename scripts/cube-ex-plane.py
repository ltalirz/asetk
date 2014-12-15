#!/usr/bin/env python
import numpy as np
import argparse
from asetk.format.cube import Cube

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
    metavar='DIRECTION',
    type=str,
    help='Plane normal. May be "x", "y" or "z"')
parser.add_argument(
    'index',
    metavar='INDEX',
    type=int,
    help='Plane index.')
parser.add_argument(
    '--plot',
    action='store_true',
    default=False,
    help='Whether to plot plane.')
parser.add_argument(
    '--replicate',
    nargs='+',
    metavar='LIST',
    type=int,
    default=None,
    help='To replicate 3x along x and 2x along y, set --replicate 3 2.')

args = parser.parse_args()

print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane, extent = c.get_plane(dir=args.dir, i=args.index, 
       return_extent=True, replica=args.replicate)

outfile = '{f}.plane{i}.dat'.format(f=args.cube, i=args.index)
print("Writing plane data to {}".format(outfile))
np.savetxt(outfile, plane)

if args.plot:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure()

    cax = plt.imshow(plane.swapaxes(0,1), extent=extent, cmap='gray')
    plt.xlabel('x [$\AA$]')
    plt.ylabel('y [$\AA$]')

    cbar = fig.colorbar(cax, format='%.2e')
    #cbar.set_label('$|\psi|^2$ $[e/a_0^2]$')
    #elif kind == 'i':
    #    cbar = fig.colorbar(cax, format='%.2f')
    #    cbar.set_label('z [$\AA$]')

    outfile = '{f}.plane{i}.png'.format(f=args.cube, i=args.index)
    print("Plotting into {}".format(outfile))
    plt.savefig(outfile, dpi=200)
