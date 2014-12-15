"""Classes for use with the Gaussian Cube file format

Provides a Cube class with reading and writing functions
"""

from __future__ import division
import numpy as np
import copy  as cp
import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as constants


class Cube(object):
    """Stores data of a cube file.
    
    Format specification according to GAUSSIAN 98:
    
    LINE   FORMAT      CONTENTS
    ===============================================================
     1     A           TITLE
     2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
     3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
     4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
     #ATOMS LINES OF ATOM COORDINATES:
     ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
     REST: 6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN
                       Y AND THEN X)
    
    FOR ORBITAL CUBE FILES, #ATOMS WILL BE < 0 AND THERE WILL BE ONE
    ADDITIONAL LINE AFTER THE FINAL ATOM GIVING THE NUMBER OF ORBITALS
    AND THEIR RESPECTIVE NUMBERS. ALSO THE ORBITAL NUMBER WILL BE
    THE FASTEST MOVING INCREMENT.
    
    ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
    """

    def __init__(self, filename=None, title=None, comment=None, origin=None, 
                 atoms=None, data=None):
        """Standard constructur, all parameters default to None."""
        self.filename = filename
        self.title = title
        self.comment = comment
        self.origin = origin
        self.atoms = atoms
        self.data = data
        self.shape_ = None   # stores shape, if grid isn't read

    @classmethod
    def from_cube(cls, cube):
        """Creates Cube from Cube object using deep copy"""
        tmp = cp.deepcopy(cube)
        return tmp

    @classmethod
    def from_file(cls, fname, read_data=False):
        """Creates Cube from cube file"""
        tmp = Cube()
        tmp.read_cube_file(fname, read_data=read_data)
        return tmp

    @property
    def cell(self):
        """Returns cell of grid.

        Unit vectors are cell[0], cell[1], cell[2].
        """
        return self.atoms.cell


    @cell.setter
    def cell(self, c):
        self.atoms.cell = c

    @property
    def shape(self):
        return self.data.shape
 
    @property
    def nx(self):
        if self.data is None:
            return self.shape_[0]
        else:
            return self.shape[0]

    @property
    def ny(self):
        if self.data is None:
            return self.shape_[1]
        else:
            return self.shape[1]

    @property
    def nz(self):
        if self.data is None:
            return self.shape_[2]
        else:
            return self.shape[2]

    @property
    def dx(self):
        return self.cell[0] / self.nx

    @property
    def dy(self):
        return self.cell[1] / self.ny

    @property
    def dz(self):
        return self.cell[2] / self.nz

    @property
    def dv(self):
        """Volume element in Angstrom^3"""
        return np.dot(self.dx, np.cross(self.dy, self.dz))

    def __str__(self):
        text  = "Spectrum containing {} spins\n".format(len(self.energylevels))
        for i in range(len(self.energylevels)):
            e = self.energylevels[i]
            s = self.spins[i]
            text += 'spin {} : {}\n'.format(s, e.__str__())
        return text

    def read_cube_file(self, fname, read_data=False, v=1):
        """Reads header and/or data of cube file
        
        """
        self.filename = fname

        f = open(fname, 'r')
        readline = f.readline

        self.title = readline()
        self.comment = readline()

        axes = [0, 1, 2]
        line = readline().split()
        natoms = int(line[0])
        b2A = constants.a0 / constants.Angstrom
        self.origin = np.array(line[1:], dtype=float) * b2A

        shape = np.empty(3)
        cell = np.empty((3, 3))
        for i in range(3):
            n, x, y, z = [float(s) for s in readline().split()]
            shape[i] = int(n)
            #if n % 2 == 1:
            #    n += 1
            cell[i] = n * np.array([x, y, z])
        self.shape_ = shape
        cell = cell * b2A

        numbers = np.empty(natoms, int)
        positions = np.empty((natoms, 3))
        for i in range(natoms):
            line = readline().split()
            numbers[i] = int(line[0])
            positions[i] = [float(s) for s in line[2:]]

        positions *= b2A 
        self.atoms = fu.Atoms(numbers=numbers, positions=positions, cell=cell)

        if read_data:
            # Note:
            # This is already ~1.7x faster than ASE's version.
            # However, parsing still dominates for reasonable disk speeds
            # (parsing time = 8x reading time on 480 MB/s SSD)
            # In order to parse quickly, use read_csv from the pandas module.

            # read() pretty much maxes out the disk read speed.
            # split() takes a considerable amount of time.
            # The conversion to float is even more expensive.
            self.data = np.array(f.read().split(), dtype=float)
            self.data = self.data.reshape(shape)
            #self.data.shape = shape
            #if axes != [0, 1, 2]:

        f.close()

    def resize(self, shape):
        """Resize dimensions of grid

        If dimensions are enlarged, additional values are filled by zeros.
        Cell vectors are updated accordingly.
        """
        shape_old = self.shape
        for dim in range(3):
            self.cell[dim] = shape[dim] * self.cell[dim] / shape_old[dim]

        tmp = np.zeros(shape)
        m = np.min([shape, self.shape], axis=0)
        tmp[:m[0], :m[1],:m[2]] = self.data[:m[0], :m[1],:m[2]] 
        self.data = tmp

        # Note: numpy's resize function always adds space
        #       at the *end* of the array.
        #self.data = np.resize(self.data, shape)

    def write_cube_file(self, fname=None):
        """Writes Cube object to file
        
        """
        if fname is None:
            fname = self.filename

        f = open(fname, 'w')

        f.write(self.title)
        f.write(self.comment)

        A2b = constants.Angstrom / constants.a0 
        o = self.origin * A2b
        f.write('{:5d}{:12.6f}{:12.6f}{:12.6f}\n' \
            .format(len(self.atoms), o[0], o[1], o[2]))

        c = self.cell * A2b
        for i in range(3): 
            n = self.shape[i] 
            d = c[i] / self.shape[i]
            f.write('{:5d}{:12.6f}{:12.6f}{:12.6f}\n'.format(n, d[0], d[1], d[2]))


        positions = self.atoms.get_positions() * A2b
        numbers = self.atoms.get_atomic_numbers() 
        for Z, (x, y, z) in zip(numbers, positions): 
            f.write('{:5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n'.format(Z, 0.0, x, y, z) ) 

        # TODO: Make 8-column format according to specs
        self.data.tofile(f, sep='\n', format='%12.6e') 
        #fmt=' %12.6e'
        #for ix in range(self.shape[0]):
        #    for iy in range(self.shape[1]):
        #        for line in range(self.shape[2] // 8 ):
        #            f.write(fmt*8.format(self.data[ix, iy, line*8 : (line+1)*8]))
        #        left = self.shape[2] % 8
        #        f.write(fmt*left.format(self.data[ix)

        f.close()

    def get_index_above_atoms(self, d, verbose=False):
        """Returns z-index of plane at z=d above topmost atom
        
        d must be given in Angstroms.
        """

        zmax = np.max(self.atoms.positions[:,2])
        zplane = zmax + d
        dz = np.linalg.norm(self.dz)

        iplane = int(round(zplane / dz))
        zplanereal = iplane * dz

        if verbose:
            print("Precise height above atoms: {} Angstroms" \
                   .format(zplanereal - zmax))
        return iplane 

    def get_plane_above_atoms(self, d, verbose=False, return_extent=None,
            replica=None, resample=None):
        """Returns plane given by z=d above topmost atom
        
        d should be given in Angstroms.
        """

        iplane =  self.get_index_above_atoms(d, verbose=verbose)
        return self.get_plane('z', iplane, return_extent=return_extent,
                              replica=replica, resample=resample)

    def get_isosurface_above_atoms(self, v, zmin=0, on_grid=False):
        """Returns z-values of isosurface

        Assumptions:
        - the tip approaches from above (i.e. along -z)
        - the values above the isosurface are smaller than below
        - the tip cannot go below zmin
        
        Parameters:
        - zmin:    minimum z-value [Angstroms] that can be reached by the tip
        - on_grid: if true, no interpolation between grid points is performed
        """

        plane = np.empty(self.shape[0:2])
        missed = 0
        dZ = np.linalg.norm(self.dz)
        nZ = self.shape[2]

        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                # argmax returns index of first occurence of maximum value
                plane[i,j] = np.argmax(self.data[i,j,::-1] > v)
                # correcting for reversing the direction
                plane[i,j] = nZ - plane[i,j] - 1 

                if plane[i,j] == nZ -1 or plane[i,j] * dZ < zmin:
                    plane[i,j] = zmin
                    missed = missed + 1
                elif on_grid:
                    plane[i,j] *= dZ
                else:
                    greater = self.data[i,j,plane[i,j]]
                    smaller = self.data[i,j,plane[i,j]+1]
                    plane[i,j] = dZ * (plane[i,j] \
                                 + (greater - v)/(greater-smaller))

        ## 1st try, just iterating through the numpy array.
        ## Turns out this is more than 40x slower than the version above.
        #for i in range(self.shape[0]):
        #    for j in range(self.shape[1]):

        #        miss = True
        #        last = self.data[i,j, self.shape[2]-1]
        #        for k in range(self.shape[2] -2, -1, -1):
        #            current = self.data[i,j,k]
        #            if current >= v and last <= v \
        #            or current <= v and last >= v:

        #                if on_grid:
        #                    plane[i,j] = k
        #                else:
        #                    # linear interpolation
        #                    plane[i,j] = dZ * (k + (last - v)/(last - current))
        #        if miss:
        #            plane[i,j] = zmin
        #            missed = missed + 1

        print("{} isovalues replaced by zmin = {}".format(missed,zmin))

        return plane

    def get_plane(self, dir, i, return_extent=None, replica=None, resample=None):
        """Returns plane normal to direction 'dir' at index 'i'
         
        * return_extent
            If True, returns [plane, extent], where extent=[x0, xmax, y0, ymax].
            Note: For plotting with matplotlib, you still need to
                  plane = plane.swapaxes(0,1)
        * replica
            replica=[3,4] will create 3x4 replicas of the original plane
        * resample
            resample=[300, 400] will resample on rectangular grid with
            300x400 points.

        """

        dvs = self.atoms.cell / self.data.shape
        ls = [ np.linalg.norm(v) for v in self.atoms.cell ]
        o = self.origin

        if dir is 'x':
            plane = self.data[i, :, :]
            dum, pdx, pdy = dvs
            pextent = [o[1], o[1]+ls[1], o[2], o[2]+ls[2]]
        elif dir is 'y':
            plane = self.data[:, i, :]
            pdy, dum, pdx = dvs
            pextent = [o[2], o[2]+ls[2], o[0], o[0]+ls[0]]
        elif dir is 'z':
            plane = self.data[:, :, i]
            dx, dy, dum = dvs
            pextent = [o[0], o[0]+ls[0], o[1], o[1]+ls[1]]
        else:
            print("Cannot recognize direction '{}'".format(dir))
            print("Direction must be 'x', 'y' or 'z'.")

        if resample:
            line = plane.flatten()
            pnx, pny = plane.shape
            pos   = [ i*pdx+j*pdy for i in range(pnx) for j in range(pny) ]
            x,y,z = zip(*pos)
         
            pextent = (np.min(x),np.max(x),np.min(y),np.max(y))
            xnew = np.linspace(pextent[0], pextent[1], resample[0])
            ynew = np.linspace(pextent[2], pextent[3], resample[1])
         
            # default interp='nn' needs mpl_toolkits.natgrid,
            # which doesn't work on some machines
            resampled = mlab.griddata(x, y, plane.flatten(), 
                                      xnew, ynew, interp='linear')
            # for some reason, I need to revert the x axis for imshow
            plane = resampled[::-1,:]

        if replica:
            plane = np.tile(plane, replica)
            pextent[1] = (pextent[1] - pextent[0]) * replica[0]
            pextent[3] = (pextent[3] - pextent[2]) * replica[1]

        if return_extent or replica or resample:
            # matplotlib will plot the 1st index along y
            # and the 2nd index along x
            #plane = plane.swapaxes(0,1)

            return [plane, pextent]
        else:
            return plane

    def set_plane(self, dir, i, plane):
        """Sets plane normal to direction 'dir' at index 'i' """

        nx, ny, nz = self.nx, self.ny, self.nz
        npx, npy = plane.shape

        if dir is 'x' and npx == ny and npy == nz:
            self.data[i, :, :] = plane
        elif dir is 'y' and npx == nz and npy == nx:
            self.data[:, i, :] = plane
        elif dir is 'z' and npx == nx and npy == ny:
            self.data[:, :, i] = plane
        else:
            print("Direction '{}' and shape '{}' are not compatible."\
                    .format(dir, plane.shape))
            print("Direction must be 'x', 'y' or 'z'.")
            return False

        return True
         
    def get_avg(self, dir):
        """Returns average value of cube file along direction 'dir'."""

        if dir is 'x':
            return np.mean(self.data, axis=0)
        elif dir is 'y':
            return np.mean(self.data, axis=1)
        elif dir is 'z':
            return np.mean(self.data, axis=2)
        else:
            print("Cannot recognize direction '{}'".format(dir))
            print("Direction must be 'x', 'y' or 'z'.")


