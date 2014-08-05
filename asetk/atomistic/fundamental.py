"""Classes of fundamental atomistic concepts

Representations of fundamental atomistic concepts,
such as energy levels, Atoms,...
Uses ASE's Atoms objects.
"""

import numpy as np
import copy  as cp
from ase.atoms import Atom, Atoms

class EnergyLevel(object):

    """An energy level"""

    def __init__(self, energy, occupation=1.0, wfn=None):
        """Set energy and occupation."""
        self.energy = energy
        self.occupation = occupation
        self.wfn = wfn

class EnergyLevels(object):

    """A list of levels with fermi energy

    Compared to the C++ class, 
    - the copy constructor has been renamed to the 'copy' method
    - the get/set methods have been stripped off
    - instead of overloading *=, the levels are exposed to the user
    - removed setFermiZero
    """

    def __init__(self, energies=None, occupations=None, wfns=None, fermi=None):
        """Fill object with levels and Fermi energy."""
        self.levels = []
        self.fermi = fermi

        # If occupations are specified, take them
        if energies is not None and occupations is not None:
            for i in range(len(energies)):
                tmp = EnergyLevel(energies[i], occupations[i])
                self.levels.append(tmp)
        # If we just have a fermi energy, create occupations
        elif energies is not None and fermi is not None:
            for i in range(len(energies)):
                if energies[i] < fermi:
                    tmp = EnergyLevel(energies[i], 1.0)
                else:
                    tmp = EnergyLevel(energies[i], 0.0)

                self.levels.append(tmp)
        # If neither fermi nor occupations are set...
        elif energies is not None:
            for i in range(len(energies)):
                tmp = EnergyLevel(energies[i], None)
                self.levels.append(tmp)

        if wfns is not None:
            if len(wfns) == len(self.levels):
                for i in range(len(wfns)):
                    self.levels[i].wfn = wfns[i]
            else:
                print("Error: Number of wave functions != number or levels")
             

    @property
    def energies(self):
        return np.array([l.energy for l in self.levels])

    @property
    def occupations(self):
        return np.array([l.occupation for l in self.levels])

    @energies.setter
    def energies(self, es):
        """Sets levels, not touching occupations."""
        if len(self.levels) != len(es):
            print('Error: Trying to set {le} energies for {le} levels.' \
                  .format(lo=len(es), le=len(self.levels)))
            return

        for i in range(len(self.levels)):
            self.levels[i].energy = es[i]

    @occupations.setter
    def occupations(self, os):
        """Sets occupations for existing levels."""

        if len(os) != len(self.levels):
            print('Error: Trying to set {lo} occupations for {le} levels.' \
                  .format(lo=len(os), le=len(self.levels)))
        else:
            for i in range(len(self.levels)):
                self.levels[i].occupation = os[i]

    def copy(self, energylevels):
        """Return a copy of energylevels."""
        self.levels = cp.copy(energylevels.levels)
        self.fermi = energylevels.fermi

    def join(self, energylevels):
        self.levels = self.levels + energylevels.levels
        if self.fermi != energylevels.fermi:
            print('Warning: Merging energy levels'
                  'with different Fermi energies.')
        self.fermi = energylevels.fermi

    def sort(self):
        self.levels.sort(key = lambda x: x.energy)

    def shift(self, de):
        self.energies += de
        self.fermi  += de

    def __iadd__(self, de):
        self.energies += de
        self.fermi  += de
        return self

    def __isub__(self, de):
        self.energies -= de
        self.fermi  -= de
        return self

    def n_occupied(self):
        """Return number of occupied levels"""

        if self.fermi:
            return sum([1 if e < self.fermi else 0 for e in self.energies])
        elif all(o is not None for o in self.occupations):
            print("Note: Counting occupied levels based on occupation number.")
            return sum([1 if o > 0 else 0 for o in self.occupations])
        else:
            print("Error: Cannot determine occupations.")

    def n_empty(self):
        """Return number of empty levels"""
        return len(levels) - self.n_occupied()

    def __str__(self):
        text  = "{} energy levels".format(len(self.levels))
        if self.fermi:
            text += ", Fermi energy {} eV".format(self.fermi)
        if all(o is not None for o in self.occupations):
            text += ", occupations specified"
        return text

    def __getitem__(self, index):
        return self.levels[index]

    def dos(self, sigma = 0.05, deltaE = 0.005, nsigma = 10):
        """
        Returns [energy, density of states].
        
        Parameters
        ----------
        sigma :  Width of Gaussian broadening [eV]
        deltaE :  spacing of energy grid [eV]
        nsigma : Gaussian distribution is cut off after nsigma*sigma
        """

        if sigma/deltaE < 10:
            print("Warning: sigma/deltaE < 10. Gaussians might not be sampled well.")
        gaussian = lambda x: 1/(np.sqrt(2*np.pi)*sigma) \
                             * np.exp(-x**2/(2*sigma**2))

        # Tabulate the Gaussian in range sigma * nsigma
        limit = sigma * nsigma
        genergies = np.r_[-limit:limit:deltaE]
        gprofile = gaussian(genergies)

        energies = self.energies
        loE=energies[0] - nsigma * sigma
        hiE=energies[-1] + nsigma * sigma
        E=np.r_[loE:hiE:deltaE]

        # Encoding the discretized energy in the array index i makes the code much faster.
        # Create dos of delta-peaks to be folded with Gaussian
        DOSdelta = np.array([0.0 for j in E])
        for e in energies:
            # In order to be able to fold with tabulated Gaussian, we have to place
            # levels *on* the grid. I.e. level spacing cannot be smaller than deltaE.
            n = int((e-loE)/deltaE)
            # Note: DOS should be calculated for unoccupied levels as well!
            #if o is not None:
            #    DOSdelta[n] += o
            #else:
            DOSdelta[n] += 1
        # Convolve with gaussian, keeping same dimension
        # Can be made even faster by using fftconvolve
        DOS = np.convolve(DOSdelta,gprofile, mode='same')

        return np.array([E, DOS])

#class Atom(object):
#
#    """Represents a single atom.
#
#    An atom with coordinates, chemical species, charge, etc.
#    """
#
#    def __init__(self, coordinates=None, number=None, charge=None):
#        self.coordinates = np.array(coordinates)
#        self.number = number
#        self.charge = charge
#
#    @property
#    def symbol(self):
#        return number2symbol(self.number)
#
#    @symbol.setter
#    def symbol(self, symbol_new):
#        self.symbol = symbol_new
#
#    def distance(self, atom):
#        d = 0
#        #print c1.coordinates
#        #print c2.coordinates
#        for c1, c2 in zip(self.coordinates, atom.coordinates):
#            d += (c1 - c2)**2
#        return np.sqrt(d)
#
#    def __str__(self):
#        text = "(x,y,z) = ({},{},{})\tN = {}\tCharge {}".format(
#                self.coordinates[0],
#                self.coordinates[1],
#                self.coordinates[2],
#                self.number,
#                self.charge)
#        return text


