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

    def __iadd__(self, b):
        if isinstance(b, float):
            self.energies += de
            self.fermi  += de
        elif isinstance(b, self.__class__):
            self.levels += b.levels
            self.sort()
            if not self.fermi == b.fermi:
                self.fermi = None
        else:
            raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'"\
                  .format(type(self), type(b)))

        return self

    def __isub__(self, de):
        self.energies -= de
        self.fermi  -= de
        return self

    def n_occupied(self, epsilon=1e-12):
        """Return number of occupied levels
        
        Levels with energy at most epsilon above Fermi  still count as occupied.
        This prevents the disregard of occupied levels, when the Fermi energy 
        is identical to the highest occupied level.
        """

        if self.fermi:
            return sum([1 if e < self.fermi + epsilon else 0 for e in self.energies])
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

    def dos(self, bmethod = 'Gaussian', bepsilon = 1e-3, FWHM = 0.1, delta_e = 0.005):
        """
        Returns [energy, density of states].
        
        Parameters
        ----------
        bmethod: Method used for broadening ('Gaussian' or 'Lorentzian')
        bepsilon: Convolution is performed with broadening function in range
            [E-Eb, E+Eb] where Eb is determined such that the integrated weight
            of the broadening function outside of [-Eb, Eb] is < bepsilon.
        FWHM: Full-width of broadening function at half-maximum [eV]
        delta_e: spacing of energy grid [eV]
        """

        import scipy.special as scsp
        # Prepare broadening functions for later convolution
        # quantile function quantile(y) = x is defined such that
        # the integral of the probability density from -\infty to x equals y.
        if bmethod == 'Gaussian':
            sigma = FWHM / np.log(8 * np.sqrt(2))
            quantile = lambda y: np.sqrt(2) * scsp.erfinv(2.0*y -1.0) * sigma
            broadening = lambda x: 1/(sigma * np.sqrt(2*np.pi)) \
                                 * np.exp( - x**2 / (2 * sigma**2) )
        elif bmethod == 'Lorentzian':
            gamma = FWHM * 0.5
            quantile = lambda y: np.tan(np.pi * (y - 0.5)) * gamma
            broadening = lambda x: 1/np.pi * gamma / (x**2 + gamma**2)
        else:
            print("Error: Broadening method \"{}\" not recognized."\
                  .format(bmethod))

        eb = -quantile(bepsilon * 0.5)
        if FWHM / delta_e < 5:
            print("Warning: FWHM / delta_e < 5. Broadening function might not be sampled well.")

        # Tabulate the Gaussian in range sigma * nsigma
        benergies = np.r_[-eb : eb : delta_e]
        bprofile = broadening(benergies)

        self.sort()
        energies = self.energies
        loE=energies[0] - eb
        hiE=energies[-1] + eb
        E=np.r_[loE : hiE : delta_e]

        # Encoding the discretized energy in the array index i makes the code much faster.
        # Create dos of delta-peaks to be folded with Gaussian
        DOSdelta = np.array([0.0 for j in E])
        for e in energies:
            # In order to be able to fold with tabulated Gaussian, we have to place
            # levels *on* the grid. I.e. level spacing cannot be smaller than deltaE.
            n = int((e-loE)/delta_e)
         
            # Note: DOS should be calculated for unoccupied levels as well!
            #if o is not None:
            #    DOSdelta[n] += o
            #else:
            DOSdelta[n] += 1
        # Convolve with gaussian, keeping same dimension
        # Can be made even faster by using fftconvolve
        DOS = np.convolve(DOSdelta,bprofile, mode='same')

        return np.array([E, DOS])


class KPoint(object):
    """Holds a k-point"""

    def __init__(self, kvector=None, energylevels=None, weight=None):
        self.kvector = np.array(kvector)
        self.energylevels = energylevels
        self.weight = weight

    def __iadd__(self, k):
        """Merging two k-points
 
        Used e.g. when merging k-points of different spins
        """
        if self.kvector != k.kvector:
            print("Warning: Adding k-points with differing k-vectors {} and {}"\
                  .format(self.kvector, k.kvector))
        self.weight += k.weight
        self.energylevels += k.energylevels 

    @property
    def nbnd(self):
        return len(self.energylevels.energies)
 
    def copy(self, kpt):
        """Performs deep copy."""
        self.energylevels = kpt.energylevels.copy()
        self.kvector = cp.copy(spectrum.kvector)
        self.weight = cp.copy(spectrum.weight)

    def __str__(self):
        e = self.energylevels
        k = self.kvector
        text = 'k = ({:6.3f}, {:6.3f}, {:6.3f})'.format(k[0], k[1], k[2])

        if self.weight:
            w = self.weight
            text += ', w = {}'.format(w)

        text += ' : {}\n'.format(e.__str__())

        return text


class Dispersion(object):
    """Holds a collection of k-points"""

    def __init__(self, kpoints=None):
         if kpoints is None:
             self.kpoints = []
         else:
             self.kpoints = kpoints

    @property
    def energylevels(self):
        s = EnergyLevels()
        for kpt in self.kpoints:
            s += kpt.energylevels
        return s

    @property
    def energies(self):
        return self.energylevels.energies

    @property
    def kvectors(self):
        s = [] 
        for kpt in self.kpoints:
            s.append(kpt.kvector)
        return s

    @property
    def weights(self):
        s = [] 
        for kpt in self.kpoints:
            s.append(kpt.weights)
        return s

    @property
    def fermi(self):
        """Returns Fermi energy."""
        fermis = [k.energylevels.fermi for k in self.kpoints]

        fermi = np.unique(fermis)

        if len( np.unique(fermis) ) != 1:
            print("There are Fermi energies {}".format(fermis))
            print("Using the mean {}".format(np.mean(fermis)))

        return np.mean(fermis)


    def merge_kpoints(self):
        kv = [0,0,0]
        levels = EnergyLevels([k.energylevels for k in self.kpoints])
        weight = 1.0
        self.kpoints = [KPoint(kv,levels,weight)]

    def __iadd__(self, s):
        """Merging two dispersions

        Used e.g. when merging dispersions belonging to different spins.
        """
        if len(self.kpoints) != len(s.kpoints):
            print("Unable to merge due to different number of kpoints")
        for i in range(len(self.kpoints)):
            self.kpoints[i] += s.kpoints[i]
        return self

    def shift(self, e):
        """Shift energylevels by energy e"""
        for k in self.kpoints:
            k.energylevels.shift(e)

    @property
    def nbnd(self):
        nbnds = [k.nbnd for k in self.kpoints]
        nbnd = np.unique(nbnds)

        if len( np.unique(nbnd) ) != 1:
            print("Warning: k-points have different numer of bands {}"\
                   .format(nbnd))
        return nbnd[0]

    @property
    def nkpt(self):
        return len(self.kpoints)


    def copy(self, dispersion):
        """Performs deep copy."""
        self.kpoints = [k.copy() for k in dispersion.kpoints]

    def __str__(self):
        text  = "Dispersion containing {} k-points\n".format(self.nkpt)
        for kpt in self.kpoints:
            text += kpt.__str__()
        return text

    def __getitem__(self, index):
        return self.kpoints[index]

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


