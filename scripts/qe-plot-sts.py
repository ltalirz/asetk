#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import asetk.format.cube as cube
import ase.io as io
import re
import argparse
import os

# Define command line parser
parser = argparse.ArgumentParser(
    description='Plots STS at given plane from QE cubefile.')
parser.add_argument('--version', action='version', version='%(prog)s 04.12.2013')
parser.add_argument(
    'cubes',
    nargs='+',
    metavar='cube files',
    help='The cube files to be sliced')
parser.add_argument(
    '--dz',
    metavar='delta z',
    default=2.5,
    type=float,
    help='The height above the topmost atom to extract the plane.')
parser.add_argument(
    '--vac',
    default=0.0,
    metavar='0.0 [eV]',
    help='The vacuum level. If specified, energies are given wrt the vacuum level.')
parser.add_argument(
    '--rep',
    default=None,
    nargs='+',
    type=int, 
    metavar='nx ny',
    help='Number of replica along x and y. If just one number is specified, it is taken for both x and y.')
parser.add_argument(
    '--samples',
    default=1000,
    metavar=int,
    help='Number of samples for resampling.')
parser.add_argument(
    '--vmax',
    default=None,
    metavar=float,
    help='If specified, maximum value of color scale for all images.')

args = parser.parse_args()

if args.rep is not None:
    if len(args.rep) == 1:
        args.rep = [ args.rep, args.rep]
    elif len(args.rep) !=2:
        print('Invalid number of replicas requested')

def get_energy(filename):
    """
    Extract energy from filename sts.-1.000.cube
    """
    e = re.search('(\-?\d+[\.\d]*?)\.cube', filename).group(1)
    return float(e) - float(args.vac)

def set_lim():
    #ring_center = [33.6,19.1]
    #ring_diameter = 38.6
    #plt.xlim([ ring_center[0]-ring_diameter/2, ring_center[0]+ring_diameter/2])
    #plt.ylim([ ring_center[1]-ring_diameter/2, ring_center[1]+ring_diameter/2])
    a=0
    #plt.ylim([0, 26])

def resample(plane, cube, rep=None, nsamples=1000):
    """Resamples data in cartesian coordinates.

    Assumptions:
    - data.shape == [self.data.shape[i] for i in axes]
    """
    nx,ny,nz    = cube.data.shape
    dx,dy,dz    = cube.atoms.cell / cube.data.shape

    if rep:
        plane = np.tile(plane, rep)
        nx *= rep[0]
        ny *= rep[1]
        
    plane = plane.flatten()
    pos   = [ i*dx+j*dy for i in range(nx) for j in range(ny) ]
    x,y,z = zip(*pos)

    extent = (np.min(x),np.max(x),np.min(y),np.max(y))
    xnew = np.linspace(extent[0], extent[1], nsamples)
    ynew = np.linspace(extent[2], extent[3], nsamples)

    resampled = mlab.griddata(x, y , plane, xnew, ynew)
    # for some reason, I need to revert the x axis for imshow
    resampled = resampled[::-1,:]

    return [resampled, extent]


files = args.cubes
fig = plt.figure(figsize=(5,5))

# cube to get data.shape
tmpcube = cube.Cube.from_file(files[0], read_data=True)
dS = tmpcube.dx[0] * tmpcube.dy[1]
print(dS)

for f in files:
    # Reading cube files is the most time consuming part of the routine.
    # Since we need only one plane out of each cube file,
    # we save it to disk for reuse.
    planefile = "{f}.z{d}".format(f=f,d=args.dz)

    plane = None
    if( os.path.isfile(planefile) ):
        plane = np.genfromtxt(planefile)
    else:
        c = cube.Cube.from_file(f, read_data=True)
        plane = c.get_plane_above_atoms(args.dz)
        # For STS at zero temperature, 
        # the occupation of the level in the calculation is irrelevant
        #plane = plane * tmp.occupation
        np.savetxt(planefile, plane)

    resampled, extent = resample(plane, cube=tmpcube, rep=args.rep, nsamples=args.samples)
    yshift = 0
    extent = [ extent[0], extent[1], extent[2] + yshift, extent[3] + yshift]

    weight = np.sum ( np.sum( plane ) ) * dS * args.rep[0] 

    print('max {m}, min {min}'.format(m=np.max(resampled), min=np.min(resampled)))

    # plotting
    plt.clf()

    cax = plt.imshow(resampled, extent=extent, vmax = float(args.vmax), cmap='gray')
    plt.xlabel('x [$\AA$]')
    plt.ylabel('y [$\AA$]')
    plt.title('E={:4.2f} eV, w = {:.1e}'.format(get_energy(f), weight))

    cbar = fig.colorbar(cax, format='%.2e')
    cbar.set_label('$\\rho(E)$')
    set_lim()
   
    #io.write('atoms.png',atoms)
    #model = plt.imread('atoms.png')
    #plt.imshow(model, extent=extent)
    
    delta=0.25
    plt.subplots_adjust(right=1-delta)
    outname='slice_{d:.2f}.png'.format(d=get_energy(f))
    plt.savefig(outname, dpi=300)
    print('Done with {f}'.format(f=f))
  
#plt.show()
   
 
    
