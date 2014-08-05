"""Classes for use with the XYZ format

Provides a Xyz class with reading and writing functions
"""

import copy  as cp
from ase.atoms import Atoms


class Xyz(object):
    """Stores data of an xyz file.
    """

    def __init__(self, filename=None, comment=None, 
                 atoms=None):
        """Standard constructor, all parameters default to None."""
        self.filename = filename
        self.comment = comment
        self.atoms = atoms

    @classmethod
    def from_xyz(cls, xyz):
        """Creates Xyz from Xyz object using deep copy"""
        tmp = cp.deepcopy(xyz)
        return tmp

    @classmethod
    def from_file(cls, fname):
        """Creates Xyz from .xyz file"""
        tmp = Xyz()
        tmp.read(fname)
        return tmp

    def read(self, fileobj, index=None):
        """ Reads a trajectory from an xyz file.

        This function returns a list of ase.atoms.Atoms objects.
        The original ase.io.xyz.read_xyz function was only capable
        of returning a single Atoms object.
        """
        if isinstance(fileobj, str):
            fileobj = open(fileobj)

        lines = fileobj.readlines()
        L1 = lines[0].split()
        if len(L1) == 1:
            del lines[:2]
            natoms = int(L1[0])
        else:
            natoms = len(lines)

        images = []
        while len(lines) >= natoms:
            positions = []
            symbols = []
            for line in lines[:natoms]:
                symbol, x, y, z = line.split()[:4]
                symbol = symbol.lower().capitalize()
                symbols.append(symbol)
                positions.append([float(x), float(y), float(z)])
            images.append(Atoms(symbols=symbols, positions=positions))
            del lines[:natoms + 2]

        if index is None:
            return images
        else:
            return images[index]

    def string(self):
        """Construct and return string for .xyz file"""
        s  = ''
        s += '  {n}\n'.format(n=len(self.atoms))
        s += self.comment + '\n'
        for atom in self.atoms:
            s += '{s:4} {x:<16.10} {y:<16.10} {z:<16.10}\n' \
                 .format(s=atom.s, x=atom.x, y=atom.y, z=atom.z)

        return s

    def write(fname=None):
        """Write content of object to .xyz file"""
        if fname is not None:
            self.fname = fname
        f=open(self.outname, 'w')
        f.write(self.string())
        f.close()

