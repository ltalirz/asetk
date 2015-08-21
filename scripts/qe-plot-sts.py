#!/usr/bin/env python
from __future__ import division
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import asetk.format.cube as cube
import asetk.format.igor as igor
import re


# Define command line parser
parser = argparse.ArgumentParser(
    description='Plots STS at given height above atoms from QE STS cubefile.')
parser.add_argument('--version', action='version', version='%(prog)s 19.05.2015')
parser.add_argument(
    '--cubes',
    nargs='+',
    metavar='FILENAME',
    help='Cube files')
#parser.add_argument(
#    '--normal',
#    nargs='+',
#    metavar='DIRECTION',
#    default='z',
#    help='Direction of the plane-normal. May be "x", "y" or "z".')
parser.add_argument(
    '--heights',
    nargs='+',
    metavar='HEIGHT',
    type=float,
    help='Height(s) of plane(s) above atoms [Angstroms].')
parser.add_argument(
    '--replicate',
    default=None,
    nargs=2,
    type=int, 
    metavar='INT',
    help='Number of replica along x and y.\
          If just one number is specified, it is taken for both x and y.')
parser.add_argument(
    '--stride',
    default=None,
    nargs=2,
    type=float, 
    metavar='INT',
    help='If specified, the data will be resampled on a cartesian grid. \
          --stride 0.5 0.5 will result in a grid twice as fine as the  \
          original grid of the cube file.')
parser.add_argument(
    '--resample',
    default=None,
    nargs=2,
    type=int, 
    metavar='INT',
    help='If specified, the data will be resampled on a cartesian grid of \
          nx x ny points.')
parser.add_argument(
    '--format',
    metavar='STRING',
    default='plain',
    help='Specifies format of output file. Can be \'plain\' (matrix of numbers)\
          or \'igor\' (igor text format of Igor Pro).'
)
parser.add_argument(
    '--plot',
    dest='plot',
    action='store_true',
    default=True,
    help='Plot data using matplotlib.')
parser.add_argument(
    '--plotrange',
    nargs=2,
    metavar='VALUE',
    default=None,
    type=float,
    help='If specified, color scale in plot will range from 1st value \
          to 2nd value.')

args = parser.parse_args()

def get_energy(filename):
    """
    Extract energy from filename sts.-1.000.cube
                              or   V_-1.000.cube
    """
    e = re.search('(\-?\d+\.?\d*?)\.cube', filename).group(1)
    #return float(e) - float(args.vac)
    return float(e)

# Make list of jobs
if args.replicate is not None:
    if len(args.replicate) == 1:
        args.replicate = [ args.replicate, args.replicate]
    elif len(args.replicate) !=2:
        print('Invalid specification of replicas. \
               Please specify --replicate <nx> <ny>.')

if args.stride is not None and args.resample is not None:
    print("Error: Please specify either --stride or --resample")
        
# Iterate over supplied cube files
for fname in args.cubes:
    print("\nReading {n} ".format(n=fname))
    c = cube.Cube.from_file(fname, read_data=True)
    dS = c.dx[0] * c.dy[1]

    if args.stride:
        s = args.stride
        resample = [ int(round(c.nx/s[0])), 
                     int(round(c.ny/s[1])) ]
    else:
        resample = args.resample

    for p in args.heights:
        planefile = None
        header = "STS at based on {}".format(fname)

        planefile = "{}.d{}".format(fname,p)
        header += ", height = {:.2f} A".format(p)
         
        plane = None
        #index = c.get_index(args.normal, p)
        plane = c.get_plane_above_atoms(p, return_object=True,
                                        replica=args.replicate, verbose=True)

        #plane = c.get_plane(args.normal, index,
        #        return_object=True, replica=args.replicate, resample=resample)

        # for details of plane object, see asetk/format/cube.py

        data = plane.data
        weight = np.sum ( np.sum( plane.data ) ) * dS
        imdata = plane.imdata
        extent = plane.extent

        if args.format == 'plain':
            datafile = planefile + '.dat'
            print("Writing {} ".format(datafile))
            np.savetxt(datafile, data, header=header)
        elif args.format == 'igor':
            igorwave = igor.Wave2d(
                    data=data,
                    xmin=extent[0],
                    xmax=extent[1],
                    xlabel='x [Angstroms]',
                    ymin=extent[2],
                    ymax=extent[3],
                    ylabel='y [Angstroms]',
            )
            datafile = planefile + '.itx'
            print("Writing {} ".format(datafile))
            igorwave.write(datafile)
        else:
            print("Error: Unknown format {}.".format(args.format))

        if args.plot:
            plotfile = planefile + '.png'
            print("Plotting into {} ".format(plotfile))
            fig = plt.figure()

            vmin = None
            vmax = None
            if args.plotrange:
                vmin = args.plotrange[0]
                vmax = args.plotrange[1]

            # when approaching from below, let smaller z be brighter
            cmap = 'gray'
            # for some reason, I need to revert the x axis for imshow
            cax = plt.imshow(imdata, extent=extent, 
                             cmap=cmap, vmin=vmin, vmax=vmax)
            plt.xlabel('x [$\AA$]')
            plt.ylabel('y [$\AA$]')

            plt.title('E={:4.2f} eV, w = {:.1e}'.format(get_energy(fname), weight))

            cbar = fig.colorbar(cax, format='%.1e')
            cbar.set_label('$\\rho(E)$ $[e/a_0^3]$')

            plt.savefig(plotfile, dpi=300, bbox_inches='tight')
            plt.close(fig)
