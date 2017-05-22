"""Classes for use with IGOR Pro

"""

import re
import numpy as np
from . import cube

class Axis(object):
    """Represents an axis of an IGOR wave"""

    def __init__(self, symbol, min, delta, unit, wavename=None):
        self.symbol = symbol
        self.min = min
        self.delta = delta
        self.unit = unit
        self.wavename = wavename

    def __str__(self):
        """Prints axis in itx format

        Note: SetScale/P expects minimum value and step-size
        """
        delta = 0 if self.delta is None else self.delta
        s = "X SetScale/P {symb} {min},{delta}, \"{unit}\", {name};\n"\
              .format(symb=self.symbol, min=self.min, delta=delta,\
                      unit=self.unit, name=self.wavename)
        return s

    def read(self, string):
        """Read axis from string

        Format: 
        X SetScale/P x 0,2.01342281879195e-11,"m", data_00381_Up;
        SetScale d 0,0,"V", data_00381_Up
        """
        match = re.search("SetScale/?P? (.) ([+-\.\de]+),([+-\.\de]+),\"(\w+)\",\s*(\w+)", string)
        self.symbol = match.group(1)
        self.min = float(match.group(2))
        self.delta = float(match.group(3))
        self.unit = match.group(4)
        self.wavename = match.group(5)



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

    def read(self, fname):
        """Read IGOR wave
        
        Should work for any dimension.
        Tested so far only for 2d wave.
        """
        f=open(fname, 'r')
        content=f.read()
        f.close()

        lines = content.split("\r")

        line = lines.pop(0)
        if not line == "IGOR":
            raise IOError("Files does not begin with 'IGOR'")

        line = lines.pop(0)
        while not re.match("WAVES",line):
            line = lines.pop(0)
        match = re.search("WAVES/N=\(([\d,]+)\)\s+(.+)",line)
        grid = match.group(1).split(',')
        grid = np.array(grid, dtype=int)
        self.name = match.group(2)

        line = lines.pop(0)
        if not line == "BEGIN":
            raise IOError("Missing 'BEGIN' statement of data block")

        # read data
        datastring = ""
        line = lines.pop(0)
        while not re.match("END",line):
            datastring += line
            line = lines.pop(0)
        data = np.array(datastring.split(), dtype=float)
        self.data = data.reshape(grid)

        # read axes
        line = lines.pop(0)
        matches = re.findall("SetScale.+?(?:;|$)", line)
        self.axes = []
        for match in matches:
            ax = Axis(None,None,None,None)
            ax.read(match)
            self.axes.append(ax)

        # the rest is discarded...
        #line = lines.pop(0)
        #print(line)

    @property
    def extent(self):
        """Returns extent for plotting"""
        grid = self.data.shape
        extent = []
        for i in range(len(grid)):
            ax = self.axes[i]
            extent.append(ax.min)
            extent.append(ax.min+ax.delta*grid[i])

        return np.array(extent)


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
        xdelta = None,
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
            x = Axis(symbol='x', min=p['xmin'], delta=p['xdelta'], unit=p['xlabel'],
                    wavename=self.name)
            self.axes = [x]

    def print_data(self):
        s = ""
        for line in self.data:
            s += "{:12.6e}\n".format(float(line))

        return s
         


class Wave2d(Wave):
    """2d Igor wave"""

    default_parameters = dict(
        xmin = 0.0,
        xdelta = None,
        xmax = None,
        xlabel = 'x',
        ymin = 0.0,
        ydelta = None,
        ymax = None,
        ylabel = 'y',
    )
 
    def __init__(self, data=None, axes=None, name=None, **kwargs):
        """Initialize 2d Igor wave

        Parameters
        ----------
        
         * data 
         * name 
         * xmin, xdelta, xlabel         
         * ymin, ydelta, ylabel         
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

            nx, ny = self.data.shape
            if p['xmax'] is None:
                p['xmax'] = p['xdelta'] * nx
            elif p['xdelta'] is None:
                p['xdelta'] = p['xmax'] / nx

            if p['ymax'] is None:
                p['ymax'] = p['ydelta'] * ny
            elif p['ydelta'] is None:
                p['ydelta'] = p['ymax'] / ny

            x = Axis(symbol='x', min=p['xmin'], delta=p['xdelta'], 
                     unit=p['xlabel'], wavename=self.name)
            y = Axis(symbol='y', min=p['ymin'], delta=p['ydelta'], 
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
            delta=c.dx,
            label="x [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='y',
            min=c.origin[1],
            delta=c.dy,
            label="y [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='z',
            min=c.origin[2],
            delta=c.dz,
            label="z [Bohr]",
            name=self.name)
            )
        axes.append(Axis(
            symbol='d',
            min=np.min(d.data),
            delta=0,
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
