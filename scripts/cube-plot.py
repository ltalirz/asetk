#!/usr/bin/env python
from __future__ import division
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import asetk.format.cube as cube
import asetk.format.qe as qe
import asetk.format.igor as igor
import sys

# Define command line parser
parser = argparse.ArgumentParser(
    description='Plots plane from Gaussian Cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 19.05.2015')
parser.add_argument(
    '--cubes',
    nargs='+',
    metavar='FILENAME',
    default=[],
    help='Cube files')
parser.add_argument(
    '--qe_cubes',
    nargs='+',
    metavar='FILENAME',
    default=[],
    help='Files in QE intermediate cube file format as written by pp.x')
parser.add_argument(
    '--normal',
    nargs='+',
    metavar='DIRECTION',
    default='z',
    help='Direction of the plane-normal. May be "x", "y" or "z".')
parser.add_argument(
    '--positions',
    nargs='+',
    metavar='HEIGHT',
    type=float,
    help='Position(s) of plane(s) along plane normal [Angstroms].')
#parser.add_argument(
#    '--from_below',
#    dest='from_below',
#    action='store_true',
#    default=False,
#    help='Approach sample from below instead from above.')
#parser.add_argument(
#    '--isovalues',
#    nargs='+',
#    metavar='VALUE',
#    type=float,
#    help='Values of the isosurface for an STM-image in constant current mode\
#          [electrons/a0^3]. 1e-7 is typically a good start.')
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
    default=(1,1),
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
    help='Plot data into png using matplotlib.')
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
if args.positions is None:
    raise ValueError("Please specify --positions")

if args.normal:
    jobs +=  zip(args.positions, [args.normal for _i in range(len(args.positions))])
    if args.normal != 'z':
        raise ValueError("Only normal z currently implemented.")
if not jobs:
    raise ValueError("No isovalues/heights specified.")

if args.replicate is not None:
    if len(args.replicate) == 1:
        args.replicate = [ args.replicate, args.replicate]
    elif len(args.replicate) !=2:
        raise ValueError('Invalid specification of replicas. \
               Please specify --replicate <nx> <ny>.')

if args.stride is not None and args.resample is not None:
    raise ValueError("Please only specify either --stride or --resample")

        
# Iterate over supplied cube files
for fname in args.cubes + args.qe_cubes:
    print("\nReading {n} ".format(n=fname))

    if fname in args.cubes:
        format = 'cube'
        c = cube.Cube.from_file(fname, read_data=True)
    elif fname in args.qe_cubes:
        format = 'qe_cube'
        tmp = qe.QECube.from_file(fname, read_data=True)
        c = tmp.to_cube()

    if args.resample:
        resample = args.resample
    elif args.stride:
        s = args.stride
        resample = [ int(round(c.nx/s[0])), 
                     int(round(c.ny/s[1])) ]

    for v,kind in jobs:
        planefile = None
        header = "STM simulation based on " + fname

        planefile = "{}.d{}{}".format(fname,kind,v)
        header += ", {} = {} [A]".format(kind,v)
        #elif kind == 'i':
        #    planefile = "{f}.iso{d}".format(f=fname,d=v)
        #    header += ", isovalue {v}, zcut {z} [A]".format(v=v, z=args.zcut)
         
        plane = None
        index = c.get_index(args.normal, v)
        plane = c.get_plane(args.normal, index,
                return_object=True, replica=args.replicate, resample=resample)
        #elif kind == 'i':
        #    plane = c.get_isosurface_above_atoms(
        #            v, zcut=args.zcut, from_below=args.from_below,
        #            return_object=True, 
        #            replica=args.replicate, resample=resample)

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
            raise ValueError("Unknown format {}.".format(args.format))

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
            cmap = 'gray'
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

            plt.savefig(plotfile, dpi=300, bbox_inches='tight')
