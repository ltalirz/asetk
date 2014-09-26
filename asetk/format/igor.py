"""Classes for use with IGOR Pro

"""

#import re
#import copy  as cp
import numpy as np
#import StringIO
#import asetk.atomistic.fundamental as fu
#import asetk.atomistic.constants as atc
import cube

class Axis(object):
    """Represents an axis of an IGOR wave"""

    def __init__(self, symbol, min, max, label, name):
        self.symbol = symbol
        self.min = min
        self.max = max
        self.label = label
        self.name = name


class Wave(object):
    """A class template for IGOR waves of generic dimension"""

    def __init__(self, data=None, name=None, axes=None):
        """Initialize IGOR wave of generic dimension"""
        self.data = data
        self.name = name
        self.axes = axes

    def __str__(self):
        """Print IGOR wave"""
        s = ""
        s += "IGOR\n"

        dimstring = "("
        for i in range(len(self.data.shape)):
            dimstring += "{}, ".format(self.data.shape[i])
        dimstring = dimstring[:-2] + ")" 

        s += "WAVES/N={}  {}\n".format(dimstring, self.name)
        s += "BEGIN\n"
        for line in self.data:
            s += "{:12.6f}\n".format(line)
        s += "END\n"
        for ax in self.axes:
            s += "X SetScale {symb} {min},{max}, \"{label}\", {name};\n"\
                  .format(symb=ax.symbol, min=ax.min, max=ax.max,\
                          label=ax.label, name=self.name)
        return s

class Wave1d(Wave):
    """1d Igor wave"""
    def __init__(self, data=None, name=None, xmin=None, xmax=None,
                 xlabel=None, ylabel=None):
        """Initialize 1d IGOR wave"""
        x = Axis(symbol = 'x', min=xmin, max=xmax, label=xlabel, name=name)
        y = Axis(symbol = 'y', min=np.min(data), max=np.max(data),
                  label=ylabel, name=name)
        super(Wave1d, self).__init__(data, name, [x,y])



class Wave2d(Wave):
    """2d Igor wave"""

    @classmethod
    def from_cube(cls, cube, dir, index, fname):
        """Creates 3d Igor Wave from Gaussian Cube file"""
        tmp = Wave3d()
        tmp.read_from_cube(fname)
        return tmp


class Wave3d(Wave):
    """3d Igor wave intended for cube files (untested)"""

    @classmethod
    def from_cube_file(cls, fname):
        """Creates 3d Igor Wave from Gaussian Cube file"""
        tmp = Wave3d()
        tmp.read_from_cube(fname)
        return tmp


    def copy(self, spectrum):
        """Performs deep copy of spectrum."""
        self.energylevels = [ el.copy() for el in spectrum.energylevels ]
        self.spins = cp.copy(spectrum.spins)

    def read_from_cube(self, fname):
        """Reads 3d Igor Wave from Gaussian Cube file"""
        c = cube.from_file(fname, read_data=True)

        self.data = c.data
        self.name = c.title

        axes = []
        axes.append(Axis(
            symbol='x',
            min=c.origin[0],
            max=c.origin[0] + c.cell[0][0],
            label="x [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='y',
            min=c.origin[1],
            max=c.origin[1] + c.cell[1][1],
            label="y [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='z',
            min=c.origin[2],
            max=c.origin[2] + c.cell[2][2],
            label="z [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='d',
            min=np.min(d.data),
            max=np.max(d.data),
            label="data [Unknown]",
            name=self.name)
            )
        self.axes = axes



#class WfnCube(cube.Cube):
#    """Gaussian cube file written by CP2K
#
#    CP2K writes the index of level and spin into the
#    comment line of the cube file
#    """
#
#    def __init__(self, title=None, comment=None, origin=None, atoms=None, 
#                 data=None, spin=None, wfn=None, energy=None, occupation=None):
#        """Standard constructor, all parameters default to None.
#        
#        energy and occupation are not stored in the cube file,
#        but can be assigned by linking the cube file with the 
#        output from the calculation.
#        """
#        super(WfnCube, self).__init__(title,comment,origin,atoms,data)
#        self.spin = spin
#        self.wfn  = wfn
#        self.energy = energy
#        self.occupation = occupation
#
#    @classmethod
#    def from_file(cls, fname, read_data=False):
#        """Creates Cube from cube file"""
#        tmp = WfnCube()
#        tmp.read_cube_file(fname, read_data=read_data)
#        return tmp
#
#    def read_cube_file(self, fname, read_data=False, v=1):
#            """Reads header and/or data of cube file"""
#            super(WfnCube, self).read_cube_file(fname, read_data, v)
#
#            # CP2K stores information on the level/spin index
#            # in the comment line
#            commentregex = 'WAVEFUNCTION\s+(\d+)\s+spin\s+(\d+)'
#            match = re.search(commentregex, self.comment)
#            self.wfn = int(match.group(1))
#            self.spin = int(match.group(2))


