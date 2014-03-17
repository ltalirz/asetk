#!/usr/bin/env python
import numpy as np
import argparse
from atk.format.cube import Cube
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def resample(plane, cube, dir, nsamples=1000):
    """Resamples data in cartesian coordinates.

    Assumptions:
    - data.shape == [self.data.shape[i] for i in axes]
    """
    nx,ny,nz    = cube.data.shape
    dx,dy,dz    = cube.atoms.cell / cube.data.shape

    if dir == 'z': 
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

    else:
        print("Not implemented yet.")

# Define command line parser
parser = argparse.ArgumentParser(
    description='Plots plane from cube file.')
parser.add_argument('--version', action='version', version='%(prog)s 27.01.2014')
parser.add_argument(
    'cube',
    metavar='',
    help='Cube file to be sliced.')
parser.add_argument(
    'dir',
    metavar='direction',
    type=str,
    help='Plane normal. May be "x", "y" or "z"')
parser.add_argument(
    'index',
    metavar='',
    type=int,
    help='Plane index.')

args = parser.parse_args()

print("Reading cube file {}".format(args.cube))
c = Cube.from_file(args.cube, read_data=True)

plane = c.get_plane(dir=args.dir, i=args.index)
resampled, extent = resample(plane, c, args.dir)

fig = plt.figure()

cax = plt.imshow(resampled, extent=extent, cmap='gray')
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

#np.savetxt(outfile, plane)
