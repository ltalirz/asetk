"""Classes for use with the Gaussian Cube file format

Provides a Cube class with reading and writing functions
"""

import numpy as np
import copy  as cp
import atomistic.fundamental as fu
import atomistic.constants as constants


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

    def __init__(self, title=None, comment=None, origin=None, 
                 cell=None, shape=None, atoms=None, data=None):
        """Standard constructur, all parameters default to None."""
        self.title = title
        self.comment = comment
        self.origin = origin
        self.cell = cell
        self.shape = shape
        self.atoms = atoms
        self.data = data

    @classmethod
    def from_file(cls, fname):
        """Creates Cube from cube file"""
        tmp = Cube()
        tmp.read_cube_file(fname)
        return tmp

    @property
    def dx(self):
        return self.cell[0] / self.shape[0]

    @property
    def dy(self):
        return self.cell[1] / self.shape[1]

    @property
    def dz(self):
        return self.cell[2] / self.shape[2]

    def copy(self, cube):
        """Performs deep copy of cube file."""
        self.title = cube.title
        self.comment = cube.comment
        self.origin = cp.copy(cube.origin)
        self.cell = cp.copy(cube.cell)
        self.shape = cp.copy(cube.shape)
        self.atoms = cp.copy(cube.atoms)
        self.data = cp.copy(cube.data)

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
        f = open(fname, 'r')
        readline = f.readline

        self.title = readline()
        self.comment = readline()

        axes = [0, 1, 2]
        line = readline().split()
        natoms = int(line[0])
        b2A = constants.a0 / constants.Angstrom
        self.origin = [b2A * float(x) for x in line[1:]]

        cell = np.empty((3, 3))
        shape = []
        for i in range(3):
            n, x, y, z = [float(s) for s in readline().split()]
            shape.append(n)
            if n % 2 == 1:
                n += 1
            cell[i] = n * b2A * np.array([x, y, z])
        self.cell = cell
        self.shape = shape

        numbers = np.empty(natoms, int)
        positions = np.empty((natoms, 3))
        for i in range(natoms):
            line = readline().split()
            numbers[i] = int(line[0])
            positions[i] = [float(s) for s in line[2:]]

        positions *= b2A 
        atoms = fu.Atoms(numbers=numbers, positions=positions, cell=cell)

        if read_data:
            # Note:
            # This is already ~1.7x faster than ASE's version.
            # However, parsing still dominates for reasonable disk speeds
            # (parsing time = 8x reading time on 480 MB/s SSD)
            # In order to parse quickly, use read_csv from the pandas module.

            # read() pretty much maxes out the disk read speed.
            # split() takes a considerable amount of time.
            # The conversion to float() is even more expensive.
            self.data = np.array(f.read().split(), dtype=float)
            self.data.reshape(shape)
            if axes != [0, 1, 2]:
                self.data = self.data.transpose(axes).copy()

