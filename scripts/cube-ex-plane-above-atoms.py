#!/usr/bin/env python
import numpy as np
import argparse
from asetk.format.cube import Cube

# Define command line parser
parser = argparse.ArgumentParser(
    description='Extracts (and plots) plane from cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 15.12.2014')
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
parser.add_argument(
    '--plot',
    action='store_true',
    default=False,
    help='Whether to plot plane.')
parser.add_argument(
    '--vmax',
    metavar='VALUE',
    default=None,
    type=float,
    help='Set maximum value of color scale for plotting.')
parser.add_argument(
    '--replicate',
    nargs='+',
    metavar='LIST',
    type=int,
    default=None,
    help='To replicate 3x along x and 2x along y, set --replicate 3 2.')

args = parser.parse_args()

# Make list of jobs
print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane, extent = c.get_plane_above_atoms(args.height, return_extent=True,
                                        replica=args.replicate)


outfile = '{f}.z{i}.dat'.format(f=args.cube, i=args.height)
print("Writing plane data to {}".format(outfile))
np.savetxt(outfile, plane)

if args.plot:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.figure()

    cax = plt.imshow(plane.swapaxes(0,1), extent=extent, cmap='gray',
                      vmax=args.vmax)
    plt.xlabel('x [$\AA$]')
    plt.ylabel('y [$\AA$]')

    cbar = fig.colorbar(cax, format='%.2e')
    #cbar.set_label('$|\psi|^2$ $[e/a_0^2]$')
    #elif kind == 'i':
    #    cbar = fig.colorbar(cax, format='%.2f')
    #    cbar.set_label('z [$\AA$]')

    outfile = '{f}.z{i}.png'.format(f=args.cube, i=args.height)
    print("Plotting into {}".format(outfile))
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
