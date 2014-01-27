"""Classes for use with Yambo

Representation of a spectrum
"""

import re
import copy  as cp
import numpy as np
import StringIO
import atk.atomistic.fundamental as fu
import atk.atomistic.constants as atc
import cube

class Dispersion:
    """A Dispersion holds the k-points belonging to one spin"""

    def __init__(self, energylevels=None, kvectors=None, weights=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.energylevels = energylevels
        self.kvectors = kvectors
        self.weights = weights

    @property
    def energies(self):
        """Returns list of energy levels of all k-points."""
        list = [el.energies for el in self.energylevels]
        return np.concatenate(list)

    @property
    def occupations(self):
        """Returns list of level occupations of all k-points."""
        os = []
        for el in self.energylevels:
            os = os + el.occupations
        return os

    def copy(self, dispersion):
        """Performs deep copy of dispersion."""
        self.energylevels = [ el.copy() for el in dispersion.energylevels ]
        self.kvectors = cp.copy(spectrum.kvectors)

    def shift(self, de):
        for levels in self.energylevels:
            levels.shift(de)

    def __str__(self):
        text  = "Dispersion containing {} k-points\n".format(len(self.energylevels))
        for i in range(len(self.energylevels)):
            e = self.energylevels[i]
            k = self.kvectors[i]
            w = self.weights[i]
            text += 'k = ({:6.3f}, {:6.3f}, {:6.3f}), w = {} : {}\n'\
                    .format(k[0], k[1], k[2], w, e.__str__())
        return text

    def __getitem__(self, index):
        return self.levels[index]

class Spectrum(object):

    """A Spectrum holds the data belonging to all spins"""

    def __init__(self, energylevels=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.dispersions = [ Dispersion(energylevels) ]

    @classmethod
    def from_output(cls, fname):
        """Creates Spectrum from Yambo output file"""
        tmp = Spectrum()
        tmp.read_from_output(fname)
        return tmp


    @property
    def energies(self):
        """Returns list of energy levels of all spins."""
        list = [disp.energies for disp in self.dispersions]
        return np.concatenate(list)

    @property
    def occupations(self):
        """Returns list of level occupations of all spins."""
        os = []
        for disp in self.dispersions:
            os = os + disp.occupations
        return os

    def copy(self, spectrum):
        """Performs deep copy of spectrum."""
        self.dispersions = [ el.copy() for el in spectrum.dispersions ]
        self.spins = cp.copy(spectrum.spins)

    def shift(self, de):
        for disp in self.dispersions:
            disp.shift(de)

    def __str__(self):
        text  = "Spectrum containing {} spins\n".format(len(self.dispersions))
        for i in range(len(self.dispersions)):
            d = self.dispersions[i]
            s = self.spins[i]
            text += 'spin {} : {}\n'.format(s+1, d.__str__())
        return text

    def __getitem__(self, index):
        return self.levels[index]

    def read_from_output(self, fname):
        s = open(fname, 'r').read()

        floatregex = '-?\d+\.\d+'
        lineregex='[^\r\n]*\r?\n'
        #blanklineregex='(?:^\s*$)'
        kptregex = 'X\* K.*?: ({f})\s*({f})\s*({f}).*?weight\s*({f}){l}(.*?)[\*\[]'\
                    .format(f=floatregex,l=lineregex)
        fermiregex='Fermi Level.*?:(\s*[\-\d\.]+)'

        self.spins=[]
        self.dispersions=[]

        # No spin for the moment, but shouldn't be too difficult to extend
        for spin in [0]:

            disp = Dispersion()
            matches=re.findall(kptregex, s, re.DOTALL)
            fermi = float(re.search(fermiregex, s).group(1))

            energylevels = []
            kvectors = []
            weights = []

            for match in matches:
                kx, ky, kz, weight, ldata  = match
                kvectors.append( np.array([kx, ky, kz], dtype=float) )
                weights.append( float(weight) )

                energies = re.findall('({f})'.format(f=floatregex), ldata, re.DOTALL)
                energies = np.array(energies)
                levels = fu.EnergyLevels(energies=energies,occupations=None, fermi=fermi)
                energylevels.append(levels)
             
            disp = Dispersion(energylevels=energylevels, kvectors=kvectors, weights=weights)

            self.dispersions.append(disp)
            self.spins.append(spin)

        

