#!/usr/bin/env python
from __future__ import division
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import asetk.format.cube as cube
import asetk.format.igor as igor
import sys

# Define command line parser
parser = argparse.ArgumentParser(
    description='Calculates Scanning Tunneling Microscopy Image \
                 from Gaussian Cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 17.12.2014')
parser.add_argument(
    '--stmcubes',
    nargs='+',
    metavar='FILENAME',
    help='Cube files containing the local density of states (s-wave tip)\
          or the appropriate matrix elements (general case).')
parser.add_argument(
    '--heights',
    nargs='+',
    metavar='HEIGHT',
    type=float,
    help='Tip-height above the topmost atom for an STM-image in constant-z\
          mode [Angstroms]. 3 Angstroms is typically reasonable.')
parser.add_argument(
    '--from_below',
    dest='from_below',
    action='store_true',
    default=False,
    help='Approach sample from below instead from above.')
parser.add_argument(
    '--isovalues',
    nargs='+',
    metavar='VALUE',
    type=float,
    help='Values of the isosurface for an STM-image in constant current mode\
          [electrons/a0^3]. 1e-7 is typically a good start.')
parser.add_argument(
    '--zcut',
    metavar='HEIGHT',
    type=float,
    default=None,
    help='Minimum z-height [Angstroms] for the tip inr constant-current mode.')
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
    '--noplot',
    dest='plot',
    action='store_false',
    default=True,
    help='Suppress plotting of resulting isosurface using matplotlib.')
parser.add_argument(
    '--plotrange',
    nargs=2,
    metavar='VALUE',
    default=None,
    type=float,
    help='If specified, color scale in plot will range from 1st value \
          to 2nd value.')

args = parser.parse_args()

# Make list of jobs
jobs = []
if args.heights:
    jobs +=  zip(args.heights, ['h' for _i in range(len(args.heights))])
if args.isovalues:
    jobs += zip(args.isovalues, ['i' for _i in range(len(args.isovalues))])
if not jobs:
    print("No isovalues/heights specified. Exiting...")
    sys.exit()    

if args.replicate is not None:
    if len(args.replicate) == 1:
        args.replicate = [ args.replicate, args.replicate]
    elif len(args.replicate) !=2:
        print('Invalid specification of replicas. \
               Please specify --replicate <nx> <ny>.')

if args.stride is not None and args.resample is not None:
    print("Error: Please specify either --stride or --resample")
        
# Iterate over supplied cube files
for fname in args.stmcubes:
    print("\nReading {n} ".format(n=fname))
    c = cube.Cube.from_file(fname, read_data=True)

    if args.stride:
        s = args.stride
        resample = [ int(round(c.nx/s[0])), 
                     int(round(c.ny/s[1])) ]
    else:
        resample = args.resample

    for v,kind in jobs:
        planefile = None
        header = "STM simulation based on " + fname
        if kind == 'h':
            planefile = "{f}.dz{d}".format(f=fname,d=v)
            header += ", z = {v} [A]".format(v=v)
        elif kind == 'i':
            planefile = "{f}.iso{d}".format(f=fname,d=v)
            header += ", isovalue {v}, zcut {z} [A]".format(v=v, z=args.zcut)
         
        plane = None
        if kind == 'h':
            plane = c.get_plane_above_atoms(v, 
                    return_object=True, from_below=args.from_below, 
                    replica=args.replicate, resample=resample)
        elif kind == 'i':
            plane = c.get_isosurface_above_atoms(
                    v, zcut=args.zcut, from_below=args.from_below,
                    return_object=True, 
                    replica=args.replicate, resample=resample)

        # for details of plane object, see asetk/format/cube.py
        data = plane.data
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
            #if kind == 'i' and args.plotrange:
            #    vmin = np.max(plane) - args.plotrange
            #    plane = plane - vmin
            #    vmin = 0

            # when approaching from below, let smaller z be brighter
            cmap = 'Greys' if args.from_below else 'gray'
            # for some reason, I need to revert the x axis for imshow
            cax = plt.imshow(imdata, extent=extent, 
                             cmap=cmap, vmin=vmin, vmax=vmax)
            plt.xlabel('x [$\AA$]')
            plt.ylabel('y [$\AA$]')

            if kind == 'h':
                cbar = fig.colorbar(cax, format='%.1e')
                cbar.set_label('$|\psi|^2$ $[e/a_0^2]$')
            elif kind == 'i':
                cbar = fig.colorbar(cax, format='%.2f')
                cbar.set_label('z [$\AA$]')

            plt.savefig(plotfile, dpi=300)
