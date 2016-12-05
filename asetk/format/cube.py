"""Classes for use with the Gaussian Cube file format

Provides a Cube class with reading and writing functions
"""

from __future__ import division
import numpy as np
import copy  as cp
import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as constants
import matplotlib.mlab as mlab

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

    dir_indices = { 
        'x': 0, 
        'y': 1, 
        'z': 2,
    }

    def __init__(self, filename=None, title=None, comment=None, origin=None, 
                 atoms=None, data=None):
        """Standard constructur, all parameters default to None."""
        self.filename = filename
        self.title = title
        self.comment = comment
        self.origin = origin
        self.atoms = atoms
        self.data = data
        self._shape = None   # stores shape, if grid data isn't read

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
        if self.data is not None:
            return self.data.shape
        elif self.shape_ is not None:
            return self.shape_
        else:
            return None
 
    @property
    def nx(self):
        return self.shape[0]

    @property
    def ny(self):
        return self.shape[1]

    @property
    def nz(self):
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

        shape = np.empty(3,dtype=int)
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

    def get_index(self, dir, d):
        """Returns index of plane at distance d along direction dir.

        d must be given in Angstroms.

        Note: Assuming orthorhombic cell.
        """
        c = self.cell

        if dir is 'x':
            dmax = c[0][0]
            step = np.linalg.norm(self.dx)
        elif dir is 'y':
            dmax = c[1][1]
            step = np.linalg.norm(self.dy)
        elif dir is 'z':
            dmax = c[2][2]
            step = np.linalg.norm(self.dz)
        else:
            raise ValueError("Did not recognize direction '{}'.".format(dir))

        if d > dmax:
            raise ValueError("Distance {} exceeds maximum distance {} \
                              along direction {}".format(d,dmax,dir))

        index = int(round(d / step))
        #d_real = index * step
        
        return index


    def get_index_above_atoms(self, d, from_below=False, verbose=False):
        """Returns z-index of plane at z=d above topmost atom
        
        d must be given in Angstroms.
        """

        if from_below:
            zmin = np.min(self.atoms.positions[:,2])
            zplane = zmmin - d
        else:
            zmax = np.max(self.atoms.positions[:,2])
            zplane = zmax + d

        dz = np.linalg.norm(self.dz)
        iplane = int(round(zplane / dz))
        zplanereal = iplane * dz

        if verbose:
            if from_below:
                delta = zmin - zplanereal
            else:
                delta = zplanereal - zmax
            print("Precise height above atoms: {} Angstroms".format(delta))

        return iplane 

    def get_plane_above_atoms(self, d, verbose=False, return_object=None,
            replica=None, resample=None, from_below=False):
        """Returns plane given by z=d above topmost atom
        
        d should be given in Angstroms.
        """

        iplane =  self.get_index_above_atoms(d, from_below, verbose=verbose)
        return self.get_plane('z', iplane, return_object=return_object,
                              replica=replica, resample=resample)

    def get_isosurface_above_atoms(self, v, from_below=False, zcut=None, 
            on_grid=False, return_object=None, replica=None, resample=None):
        """Returns z-values of isosurface

        Assumptions:
        - the values above the isosurface are smaller than below
        - the tip cannot go below zcut
        
        Parameters:
        - from_below: tip approaches from below instead of from above.
        - zcut:       minimum z-value [Angstroms] that can be reached by the tip
                      (maximum z-value for approach from below)
        - on_grid:    if true, no interpolation between grid points is performed
        """

        plane = np.empty(self.shape[0:2])
        missed = 0

        dz = np.linalg.norm(self.dz)
        nz = self.nz

        if zcut is None:
            zcut = dz*nz if from_below else 0.0

        if from_below:
            zmax = zcut
        else:
            # the following is written for an "approach from below"
            self.data = self.data[:,:,::-1]
            zmax = nz*dz - zcut

        # this is written for an "approach from below"
        for i in range(self.nx):
            for j in range(self.ny):
                # argmax returns index of first occurence of maximum value
                itmp = np.argmax(self.data[i,j,:] > v)

                if itmp == 0 or itmp * dz > zmax:
                    plane[i,j] = zmax
                    missed = missed + 1
                elif on_grid:
                    plane[i,j] = itmp * dz
                else:
                    greater = self.data[i,j,itmp]
                    smaller = self.data[i,j,itmp-1]
                    plane[i,j] = dz * (itmp - (greater - v)/(greater-smaller))

        # revert back to original data set
        if not from_below:
            self.data = self.data[:,:,::-1]
            plane = dz*nz - plane

            ##v2
            ## argmax returns index of first occurence of maximum value
            #plane[i,j] = np.argmax(self.data[i,j,::-1] > v)
            ## correcting for reversing the direction
            #plane[i,j] = nz - plane[i,j] - 1 

            #if plane[i,j] == nz -1 or plane[i,j] * dz < zmin:
            #    plane[i,j] = zmin
            #    missed = missed + 1
            #elif on_grid:
            #    plane[i,j] *= dz
            #else:
            #    greater = self.data[i,j,plane[i,j]]
            #    smaller = self.data[i,j,plane[i,j]+1]
            #    plane[i,j] = dz * (plane[i,j] \
            #                 + (greater - v)/(greater-smaller))

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

        print("{} z-values replaced by zcut = {}".format(missed,zcut))

        pextent, pdx, pdy = self.get_plane_extent('z', return_vectors=True)

        plane = Plane(data=plane, origin=self.origin, dx=pdx, dy=pdy)

        if replica:
            plane.replicate(replica)

        if resample:
            plane.resample(resample)

        if return_object:
            # matplotlib will plot the 1st index along y
            # and the 2nd index along x
            #plane = plane.swapaxes(0,1)

            return plane
        else:
            return plane.data

    def get_plane(self, dir, i, return_object=False, replica=None, resample=None):
        """Returns plane normal to direction 'dir' at index 'i'
         
        * return_object
            If True, returns Plane object (which knows about its extent).
        * replica
            replica=[3,4] will create 3x4 replicas of the original plane
        * resample
            resample=[300, 400] will resample on rectangular grid with
            300x400 points.

        """

        shape = self.data.shape
        dvs = self.atoms.cell / shape
        ls = [ np.linalg.norm(v) for v in self.atoms.cell ]
        o = self.origin

        if dir is 'x' and i < shape[0]:
            plane = self.data[i, :, :]
        elif dir is 'y' and i < shape[1]:
            plane = self.data[:, i, :]
        elif dir is 'z' and i < shape[2]:
            plane = self.data[:, :, i]
        else:
            msg  = "Direction {} not recognized or index {} out of bounds"\
                    .format(dir,i)
            msg += "\nDirection must be 'x', 'y' or 'z'."
            raise ValueError(msg)

        pextent, pdx, pdy = self.get_plane_extent(dir, return_vectors=True)

        plane = Plane(data=plane, origin=o, dx=pdx, dy=pdy)


        if replica:
            plane.replicate(replica)

        if resample:
            plane.resample(resample)

        if return_object:
            return plane
        else:
            return plane.data


    def get_plane_extent(self, dir, return_vectors=False):
        """Returns extent of plane.

        Useful for plotting with matplotlib.

         * return_vectors : If True, returns [pextent, pdx, pdy]
             where pdx, pdy are the vectors of the plane grid.
        """

        dvs = self.atoms.cell / self.data.shape
        ls = [ np.linalg.norm(v) for v in self.atoms.cell ]
        o = self.origin

        if dir is 'x':
            dum, pdx, pdy = dvs
            pextent = [o[1], o[1]+ls[1], o[2], o[2]+ls[2]]
        elif dir is 'y':
            pdy, dum, pdx = dvs
            pextent = [o[2], o[2]+ls[2], o[0], o[0]+ls[0]]
        elif dir is 'z':
            pdx, pdy, dum = dvs
            pextent = [o[0], o[0]+ls[0], o[1], o[1]+ls[1]]
        else:
            print("Cannot recognize direction '{}'".format(dir))
            print("Expected 'x', 'y' or 'z'.")

        if return_vectors:
            return [pextent, pdx, pdy]
        else:
            return pextent
        


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

    def __iadd__(self, c):
        """Adds grid values of two cube files together"""
        if self.data.shape != c.data.shape:
            raise ValueError("Shape of cube files do not agree")
        self.data += c.data

        return self

    def roll(self, dir, shift=None, distance=None):
        """Rolls cube along direction dir

        Positions of atoms are updated accordingly.

        parameters
        ----------
        dir: string
          'x', 'y' or 'z'
        shift: int
           how many indices to shift data along dir
        dist: float
           (alternative) what distance to shift along dir
        
        """
        dir_index = self.dir_indices[dir]
        step = np.linalg.norm(self.cell[dir_index] / self.data.shape[dir_index])


        if shift:
            dist_exact = shift * step
        elif distance:
            shift = int(distance / step)
            dist_exact = shift * step
        else:
            raise IOError("Please provide either shift (integer) or distance (float)")

        print("Rolling cube file by {:.3f} Angstroms along {}."\
                .format(dist_exact,dir))
        self.data = np.roll(self.data, shift=shift, axis=dir_index)

        v = np.zeros(3)
        v[dir_index] = dist_exact
        self.atoms.translate(v)


class Plane(object):
    """Stores a plane of a cube file.
    
    """

    def __init__(self, data=None, origin=None, dx=None, dy=None, extent=None):
        """Standard constructur, all parameters default to None."""
        self.data = data
        self.origin = origin

        if extent != None and (dx != None or dy != None):
            print("Error: Please specify either extent or dx, dy")
        elif extent != None:
            self.dx = [(extent[1]-extent[0])/self.nx, 0]
            self.dy = [0, (extent[3]-extent[2])/self.ny]
        else:
            self.dx = dx
            self.dy = dy

    @property
    def nx(self):
        return self.data.shape[0]

    @property
    def ny(self):
        return self.data.shape[1]

    @property
    def imdata(self):
        """Returns data in proper form for plotting with matplotlib
        
        The x-y plane of a cube file is a matrix of the form
        
        x1y1, x1y2, ...
        x2y1, x2y2, ...
        ...

        For plotting with matplotlib, we want the form

        ...
        x1y2, x2y2, ...
        x1y1, x2y1, ...

        """
        return  (self.data.swapaxes(0,1))[::-1,:]

    @property
    def extent(self):
        """Returns extent of plane.

        Useful for plotting with matplotlib.
        """
        o = self.origin
        dx = np.linalg.norm(self.dx)
        dy = np.linalg.norm(self.dy)

        extent = [ o[0], o[0]+dx*self.nx, o[1], o[1]+dy*self.ny ]

        return extent


    def replicate(self, replica):
        self.data = np.tile(self.data, replica)
        e = self.extent

        self.extent[0] = (e[1] - e[0]) * replica[0]
        self.extent[3] = (e[3] - e[2]) * replica[1]


    def resample(self, npoints):

        pdx = self.dx
        pdy = self.dy
        pnx = self.nx
        pny = self.ny
        o = self.origin

        line = self.data.flatten()
        pos   = [ o + i*pdx + j*pdy for i in range(pnx) for j in range(pny) ]
        x,y,z = zip(*pos)
     
        extent = [np.min(x),np.max(x),np.min(y),np.max(y)]
        xnew = np.linspace(extent[0], extent[1], npoints[0])
        ynew = np.linspace(extent[2], extent[3], npoints[1])
        dx = (extent[1]-extent[0]) / npoints[0]
        dy = (extent[2]-extent[3]) / npoints[1]
        self.dx = [dx, 0]
        self.dy = [0, dy]
     
        # default interp='nn' needs mpl_toolkits.natgrid,
        # which doesn't work on some machines
        resampled = mlab.griddata(x, y, line,
                                  xnew, ynew, interp='linear')
        # translating the output from mlab to floats
        resampled = np.array(resampled, dtype=float)

        # in order to go back to the original data layout
        # (1st index x, 2nd index y), we need to swap
        self.data = resampled.swapaxes(0,1)

    
