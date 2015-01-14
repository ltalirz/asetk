"""Classes for use with IGOR Pro

"""

#import re
#import copy  as cp
import numpy as np
#import StringIO
#import asetk.atomistic.fundamental as fu
#import asetk.atomistic.constants as atc
from . import cube

class Axis(object):
    """Represents an axis of an IGOR wave"""

    def __init__(self, symbol, min, max, unit, wavename=None):
        self.symbol = symbol
        self.min = min
        self.max = max
        self.unit = unit
        self.wavename = wavename

    def __str__(self):
        max = 0 if self.max is None else self.max
        s = "X SetScale {symb} {min},{max}, \"{unit}\", {name};\n"\
              .format(symb=self.symbol, min=self.min, max=max,\
                      unit=self.unit, name=self.wavename)
        return s



class Wave(object):
    """A class template for IGOR waves of generic dimension"""

    def __init__(self, data, axes, name=None):
        """Initialize IGOR wave of generic dimension"""
        self.data = data
        self.name = "PYTHON_IMPORT" if name is None else name
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
        s += self.print_data()
        s += "END\n"
        for ax in self.axes:
            s += str(ax)
        return s

    def print_data(self):
        """Determines how to print the data block.
        
        To be implemented by subclasses."""

    def write(self, fname):
        f=open(fname, 'w')
        f.write(str(self))
        f.close()


class Wave1d(Wave):
    """1d Igor wave"""
    
    default_parameters = dict(
        xmin = 0.0,
        xmax = None,
        xlabel = 'x',
        ylabel = 'y',
    )

    def __init__(self, data=None, axes=None, name="1d", **kwargs):
        """Initialize 1d IGOR wave"""
        super(Wave1d, self).__init__(data, axes, name) 

        self.parameters = self.default_parameters
        for key, value in kwargs.items():
            if key in self.parameters:
                self.parameters[key] = value
            else:
                raise KeyError("Unknown parameter {}".format(key))

        if axes is None:
            p=self.parameters
            x = Axis(symbol='x', min=p['xmin'], max=p['xmax'], label=p['xlabel'])
            self.axes = [x]

    def print_data(self):
        s = ""
        for line in self.data:
            s += "{:12.6e}\n".format(float(line))
         


class Wave2d(Wave):
    """2d Igor wave"""

    default_parameters = dict(
        xmin = 0.0,
        xmax = 1.0,
        xlabel = 'x',
        ymin = 0.0,
        ymax = 1.0,
        ylabel = 'y',
    )
 
    def __init__(self, data=None, axes=None, name=None, **kwargs):
        """Initialize 2d Igor wave

        Parameters
        ----------
        
         * data 
         * name 
         * xmin, xmax, xlabel         
         * ymin, ymax, ylabel         
        """
        super(Wave2d, self).__init__(data, axes=axes, name=name)

        self.parameters = self.default_parameters
        for key, value in kwargs.items():
            if key in self.parameters:
                self.parameters[key] = value
            else:
                raise KeyError("Unknown parameter {}".format(key))

        if axes is None:
            p=self.parameters
            x = Axis(symbol='x', min=p['xmin'], max=p['xmax'], 
                     unit=p['xlabel'], wavename=self.name)
            y = Axis(symbol='y', min=p['ymin'], max=p['ymax'], 
                     unit=p['ylabel'], wavename=self.name)
            self.axes = [x,y]


    def print_data(self):
        """Determines how to print the data block"""
        s = ""
        for line in self.data:
            for x in line:
                s += "{:12.6e} ".format(x)
            s += "\n"

        return s
         

    @classmethod
    def from_cube(cls, cube, dir, index, fname):
        """Creates 2d Igor Wave from Gaussian Cube file
        
        Parameters
        ----------
         * cube : format.cube object containing cube file
         * dir  : 'x', 'y' or 'z'
         * index: index of plane to be taken
         """
        tmp = Wave3d()
        tmp.read_from_cube(fname)
        return tmp

    def read_from_cube(self, cube, dir, index, fname=None):
        # To be implemented
        dir
        #super(Wave2d, self).__init__(
        #        data=cube.get_plane(dir, index),
        #        name=name,
        #        axes=)
        #        comment=comment,
        #        t,origin,atoms,data)
        


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


