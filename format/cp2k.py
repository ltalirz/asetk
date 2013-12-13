"""Classes for use with CP2K

Representation of a spectrum
"""

import re
import copy  as cp
import numpy as np
import StringIO
import atomistic.fundamental as fu
import atomistic.constants as atc
import cube

class Spectrum(object):

    """An collection of energy levels, grouped by spin"""

    def __init__(self, energylevels=None, spins=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.energylevels = energylevels
        self.spins = spins

    @classmethod
    def from_mo(cls, fname):
        """Creates Spectrum from cube file"""
        tmp = Spectrum()
        tmp.read_from_mo(fname)
        return tmp


    @property
    def energies(self):
        """Returns list of energy levels of all spins."""
        ens = []
        for el in self.energylevels:
            ens = ens + el.energies
        return ens

    @property
    def occupations(self):
        """Returns list of level occupations of all spins."""
        os = []
        for el in self.energylevels:
            os = os + el.occupations
        return os

    def copy(self, spectrum):
        """Performs deep copy of spectrum."""
        self.energylevels = [ el.copy() for el in spectrum.energylevels ]
        self.spins = cp.copy(spectrum.spins)

    def shift(self, de):
        self.energies -= de
        self.fermi  -= de

    def __str__(self):
        text  = "Spectrum containing {} spins\n".format(len(self.energylevels))
        for i in range(len(self.energylevels)):
            e = self.energylevels[i]
            s = self.spins[i]
            text += 'spin {} : {}\n'.format(s, e.__str__())
        return text

    def __getitem__(self, index):
        return self.levels[index]

    def read_from_mo(self, fname):
        s = open(fname, 'r').read()

        lineregex='[^\r\n]*\r?\n'
        fermiregex='.*?Fermi energy:(\s*[\-\d\.]+)'
        matches=re.findall(
                'EIGENVALUES({l}){l}{l}([\-\.\s\d]*){f}' \
                    .format(l=lineregex, f=fermiregex), 
                s, re.DOTALL)

        self.spins=[]
        self.energylevels=[]
        spin=0
        for match in matches:
            # we take only the last ones, which do not contain 'SCF'
            if not re.search('SCF', match[0]):
                data = np.genfromtxt(StringIO.StringIO(match[1]),
                             dtype=[int,float,float])
                i,E,occ = zip(*data)
                fermi = float(match[2]) * atc.Ha / atc.eV
                E     = np.array(E)     * atc.Ha / atc.eV

                levels = fu.EnergyLevels(energies=E,occupations=occ, fermi=fermi)
                self.energylevels.append(levels)
                self.spins.append(spin)

                spin = spin + 1

        

class WfnCube(cube.Cube):
    """Gaussian cube file written by CP2K

    CP2K writes the index of level and spin into the
    comment line of the cube file
    """


    def __init__(self, title=None, comment=None, origin=None,
                 cell=None, shape=None, atoms=None, data=None, 
                 spin=None, wfn=None, energy=None):
        """Standard constructur, all parameters default to None."""
        cube.Cube.__init__(self, title,comment,origin,cell,shape,atoms,data)
        self.spin = spin
        self.wfn  = wfn
        self.energy = energy

    @classmethod
    def from_file(cls, fname):
        """Creates Cube from cube file"""
        tmp = WfnCube()
        tmp.read_cube_file(fname)
        return tmp

    def read_cube_file(self, fname, read_data=False, v=1):
            """Reads header and/or data of cube file"""
            cube.Cube.read_cube_file(self, fname, read_data, v)

            # CP2K stores information on the level/spin index
            # in the comment line
            commentregex = 'WAVEFUNCTION\s+(\d+)\s+spin\s+(\d+)'
            match = re.search(commentregex, self.comment)
            self.wfn = int(match.group(1))
            self.spin = int(match.group(2))



