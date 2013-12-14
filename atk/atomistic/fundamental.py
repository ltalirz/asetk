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
            for i in xrange(len(energies)):
                tmp = EnergyLevel(energies[i], occupations[i])
                self.levels.append(tmp)
        # If we just have a fermi energy, create occupations
        elif energies is not None and fermi is not None:
            for i in xrange(len(energies)):
                if energies[i] < fermi:
                    tmp = EnergyLevel(energies[i], 1.0)
                else:
                    tmp = EnergyLevel(energies[i], 0.0)

                self.levels.append(tmp)
        # If neither fermi nor occupations are set...
        elif energies is not None:
            for i in xrange(len(energies)):
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
        """Sets levels, not specifying occupations."""
        self.levels = [ EnergyLevels(es[i], None) for i in range(len(ens)) ]

    @occupations.setter
    def occupations(self, os):
        """Sets occupations for existing levels."""

        if len(os) != len(self.levels):
            print('Error: Trying to set {lo} occupations for {le} levels.' \
                  .format(lo=len(os), le=len(self.levels)))
        else:
            for i in xrange(len(self.levels)):
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
        self.energies -= de
        self.fermi  -= de

    def n_occupied(self):
        """Return number of occupied levels"""

        if fermi:
            return sum([1 if e < fermi else 0 for e in energies])
        elif all(o is not None for o in self.occupations):
            print("Note: Counting occupied levels based on occupation number.")
            return sum([1 if o > 0 else 0 for o in occupations])
        else:
            print("Error: Cannot determine occupations.")

    def n_unoccopied(self):
        """Return number of unoccupied levels"""
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


