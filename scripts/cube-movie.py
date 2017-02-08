#!/usr/bin/env python
from __future__ import division, print_function
import numpy as np
import argparse
import re
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ase.data import covalent_radii
from ase.data.colors import cpk_colors

import asetk.format.cube as cube
import asetk.format.qe as qe

# Define command line parser
parser = argparse.ArgumentParser(
    description='Plots movie from Gaussian Cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 31.01.2017')
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
    '--sts_cubes',
    nargs='+',
    metavar='FILENAME',
    default=[],
    help='STS cube files')
parser.add_argument(
    '--normal',
    nargs='+',
    metavar='DIRECTION',
    default='z',
    help='Direction along which to make the movie. May be "x", "y" or "z".')
#parser.add_argument(
#    '--replicate',
#    default=None,
#    nargs=2,
#    type=int, 
#    metavar='INT',
#    help='Number of replica along x and y.\
#          If just one number is specified, it is taken for both x and y.')
#parser.add_argument(
#    '--stride',
#    default=(1,1),
#    nargs=2,
#    type=float, 
#    metavar='INT',
#    help='If specified, the data will be resampled on a cartesian grid. \
#          --stride 0.5 0.5 will result in a grid twice as fine as the  \
#          original grid of the cube file.')
#parser.add_argument(
#    '--resample',
#    default=None,
#    nargs=2,
#    type=int, 
#    metavar='INT',
#    help='If specified, the data will be resampled on a cartesian grid of \
#          nx x ny points.')
parser.add_argument(
    '--format',
    metavar='STRING',
    default='mp4',
    help='Specifies format of output. Can be \'png\' (collection of pngs) or \'mp4\' (movie).'
)
parser.add_argument(
    '--plotrange',
    nargs=2,
    metavar='VALUE',
    default=None,
    type=float,
    help='If specified, color scale in plot will range from 1st value \
          to 2nd value.')
parser.add_argument(
    '--ffmpeg_path',
    metavar='FILEPATH',
    default='/opt/local/bin/ffmpeg',
    help='Specify path to ffmpeg (required for mp4 movies).')
parser.add_argument(
    '--atom_sketch',
    metavar='STRING',
    default='False',
    help='If True overlay of atomic sketch'
)
parser.add_argument(
    '--atom_zcut',
    metavar='VALUE',
    default=100.0,
    type=float,
    help='cuts atoms below zcut in the atomic sketch'
)


args = parser.parse_args()

if args.format not in ['png','mp4']:
    raise ValueError("Only png and mp4 format currently supported.")

plt.rcParams['animation.ffmpeg_path'] = args.ffmpeg_path


dir_index = cube.Cube.dir_indices[args.normal]

# Iterate over supplied cube files
for fname in args.cubes + args.qe_cubes + args.sts_cubes:
    print("\nReading {n} ".format(n=fname))

    if fname in args.cubes:
        cformat = 'cube'
        c = cube.Cube.from_file(fname, read_data=True)
    elif fname in args.sts_cubes:
        cformat = 'sts_cube'
        c = cube.STSCube.from_file(fname, read_data=True)
    elif fname in args.qe_cubes:
        cformat = 'qe_cube'
        tmp = qe.QECube.from_file(fname, read_data=True)
        c = tmp.to_cube()


    name = os.path.splitext(fname)[0]

    vmax = np.amax(c.data)
    vmin = np.amin(c.data)

    fig, (ax) = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(6,5))
    plt.subplots_adjust(left=0.10)
    ax.set_xlabel('x [$\AA$]')
    ax.set_ylabel('y [$\AA$]')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ax.set_title("{}".format(c.title))

    if cformat == 'sts_cube':
        label_text = "U = {:.3f} V"
    else:
        label_text = "z = {:.3f}"

    z0 = c.origin[dir_index]
    nz = c.data.shape[dir_index]
    dz = np.linalg.norm(c.cell[dir_index]) / nz

######plot atomic sketch

    if args.atom_sketch=='True':
        if args.atom_zcut < 100 :
            xsketch=[xx[0]  for xx in c.atoms.positions if xx[2] > args.atom_zcut]
            ysketch=[xx[1]  for xx in c.atoms.positions if xx[2] > args.atom_zcut]
            nsketch=[c.atoms.numbers[i]   for i in range(len(c.atoms.positions)) if c.atoms.positions[i][2] > args.atom_zcut]
            sketch=zip(nsketch,xsketch,ysketch)

        circ=[]
        for a in sketch:
            circ.append(plt.Circle((a[1],a[2] ), covalent_radii[a[0]], color=cpk_colors[a[0]],fill=False,clip_on=True))
        for thecirc in circ:
            ax.add_artist(thecirc)
######end plot atomic sketch

    ims = []
    for i in range(nz):
        p = c.get_plane(args.normal,i,return_object=True)
        im = ax.imshow(p.imdata, norm=plt.Normalize(vmin,vmax), 
            extent=p.extent, cmap=matplotlib.cm.bwr)
        label = ax.text(0.8, 0.9,label_text.format(z0+i*dz),
             horizontalalignment='center', verticalalignment='center',
             transform = ax.transAxes)

        if i == 0:
            plt.colorbar(im, cax=cax)

        if args.format=='png':
             outname = "{}_{}_{:04d}.png".format(name,args.normal,i)
             print("Saving {}".format(outname), end='\r')
             sys.stdout.flush()
             plt.savefig(outname, dpi=200)
             label.remove()
             im.remove()
        elif args.format=='mp4':
             ims.append( (im, label,) )


    if args.format =='mp4':
        im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
                                           blit=True)
        mfile = "{}_{}.mp4".format(name,args.normal)
        print("Making movie {}".format(mfile))
        writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        #writer = animation.MencoderWriter(fps=15, bitrate=1800)
        im_ani.save(mfile, writer=writer)

