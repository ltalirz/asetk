"""Classes for use with CP2K

Representation of a spectrum
"""

import re
import copy as cp
import numpy as np
import io
import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as atc
from . import cube


class Spectrum(object):

    """An collection of energy levels, grouped by spin"""

    def __init__(self, energylevels=None, spins=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.energylevels = energylevels
        self.spins = spins

    @classmethod
    def from_mo(cls, fname):
        """Creates Spectrum from list of molecular occupations"""
        tmp = Spectrum()
        tmp.read_from_mo(fname)
        return tmp

    @classmethod
    def from_output(cls, fname):
        """Creates Spectrum from CP2K output"""
        tmp = Spectrum()
        tmp.read_from_output(fname)
        return tmp

    @classmethod
    def from_pdos(cls, fname):
        """Creates Spectrum from projected density of states"""
        tmp = Spectrum()
        tmp.read_from_pdos(fname)
        return tmp

    @property
    def energies(self):
        """Returns list of energy levels of all spins."""
        list = [el.energies for el in self.energylevels]
        return np.concatenate(list)

    @property
    def occupations(self):
        """Returns list of level occupations of all spins."""
        os = np.array([])
        for el in self.energylevels:
            os = np.concatenate((os, el.occupations))
        return os

    @property
    def fermi(self):
        """Returns Fermi energy."""
        fermis = [el.fermi for el in self.energylevels]

        if len(np.unique(fermis)) != 1:
            print(("There are Fermi energies {}".format(fermis)))
            print(("Using the mean {}".format(np.mean(fermis))))

        return np.mean(fermis)

    def copy(self, spectrum):
        """Performs deep copy of spectrum."""
        self.energylevels = [el.copy() for el in spectrum.energylevels]
        self.spins = cp.copy(spectrum.spins)

    def n_occupied(self):
        """Return number of occupied levels"""
        noccs = [el.n_occupied() for el in self.energylevels]
        if len(np.unique(noccs)) != 1:
            print("The number of occupied orbitals {} varies between spins"
                  .format(noccs))
            print("Using the mean {}".format(np.mean(noccs)))

        return int(np.mean(noccs))

    def n_empty(self):
        """Return number of empty levels"""
        nempts = [el.n_empty() for el in self.energylevels]
        if len(np.unique(nempts)) != 1:
            print("The number of empty orbitals {} varies between spins"
                  .format(nempts))
            print("Using the mean {}".format(np.mean(nempts)))

        return int(np.mean(nempts))

    def shift(self, de):
        for levels in self.energylevels:
            levels.shift(de)

    def __iadd__(self, de):
        for levels in self.energylevels:
            levels.shift(de)
        return self

    def __isub__(self, de):
        for levels in self.energylevels:
            levels.shift(-de)
        return self

    def dos(self, bmethod='Gaussian', bepsilon=1e-3, FWHM=0.1, delta_e=0.005):
        """
        Returns [energy, density of states].

        For documentation of parameters, see 
          atk.atomistic.fundamental.EnergyLevels.dos
        """

        nspin = len(self.energylevels)
        if nspin > 1:
            print("Summing contributions from {} spins"
                  .format(len(self.energylevels)))

            fermis = [el.fermi for el in self.energylevels]
            if len(np.unique(fermis)) != 1:
                print("Warning: Fermi energies {} differ".format(fermis))

            DOSes = []
            for el in self.energylevels:
                E, DOS = el.dos(bmethod, bepsilon, FWHM, delta_e)
                DOSes.append(DOS)

            return [E, np.sum(DOSes, axis=0)]

        elif nspin == 1:
            return self.energylevels[0].dos(bmethod, bepsilon, FWHM, delta_e)

        else:
            print("Error: DOS requested, but no states found.")
            return

    def __str__(self):
        text = "Spectrum containing {} spins\n".format(len(self.energylevels))
        for i in range(len(self.energylevels)):
            e = self.energylevels[i]
            s = self.spins[i]
            text += 'spin {} : {}\n'.format(s+1, e.__str__())
        return text

    def __getitem__(self, index):
        return self.levels[index]

    def read_from_mo(self, fname):
        """Reads Spectrum from list of molecular occupations"""
        s = open(fname, 'r').read()

        lineregex = '[^\r\n]*\r?\n'
        fermiregex = 'Fermi energy:(\s*[\-\d\.]+)'
        mainregex = 'MO EIGENVALUES AND MO OCCUPATION NUMBERS(?!, AND)({l}){l}{l}([\-\.\s\d]*).*?{f}' \
            .format(l=lineregex, f=fermiregex)
        spinregex = '(?:ALPHA|BETA)\s*'

        if re.search(spinregex, s):
            matches = re.findall(spinregex+mainregex, s, re.DOTALL)
            last_matches = [matches[-2], matches[-1]]
        else:
            matches = re.findall(mainregex, s, re.DOTALL)
            last_matches = [matches[-1]]

        # first, we sort out all but the last one
        #last_matches = []
        # print(matches[0])
        # for match in reversed(matches):
        #    #print(match)
        #    # we take only the ones, which do not contain 'SCF'
        #    if re.search('SCF', match[0]):
        #        continue
        #    elif re.search('ALPHA', match[0]) \
        #      or re.search('BETA', match[0]):
        #          if len(last_matches) == 0:
        #              last_matches.append(match)
        #              continue
        #          else:
        #              last_matches.append(match)
        #              break
        #    else:
        #        last_matches.append(match)
        #        break

        self.spins = []
        self.energylevels = []
        spin = 0
        for match in last_matches:
            # we take only the last ones, which do not contain 'SCF'
            # if not re.search('SCF', match[0]):

            # python2
            # data = np.genfromtxt(StringIO.StringIO(match[2]),
            #             dtype=[int,float,float])
            data = np.genfromtxt(io.BytesIO(match[1].encode()),
                                 dtype=[int, float, float])
            i, E, occ = list(zip(*data))
            fermi = float(match[2]) * atc.Ha / atc.eV
            E = np.array(E) * atc.Ha / atc.eV

            levels = fu.EnergyLevels(energies=E, occupations=occ, fermi=fermi)
            self.energylevels.append(levels)
            self.spins.append(spin)

            spin = spin + 1

    def read_from_output(self, fname):
        """Reads Spectrum from CP2K output"""
        s = open(fname, 'r').read()

        # Format in CP2K output
        if not re.search('Eigenvalues of the occupied subspace', s):
            raise ValueError("Unable to parse CP2K output file")

        # We are interested only in the last run (result of calculation)
        rundatas = re.findall(
            'Eigenvalues of the occupied subspace .*?FORCE_EVAL', s, re.DOTALL)
        # rundatas = re.findall(
        #        'Eigenvalues of the occupied subspace .*?HOMO', s, re.DOTALL)
        rundata = rundatas[-1]

        # each spin may have occupied & unoccupied levels + fermi energy
        occupieddatas = re.findall(
            'Eigenvalues of the occupied subspace spin([\.\-\d\s]*)',
            rundata, re.DOTALL)
        #unoccupieddatas = re.findall('Lowest eigenvalues of the unoccupied subspace spin.*?iterations([\.\-\d\s]*)', rundata, re.DOTALL)
        unoccupieddatas = re.findall(
            'Lowest (?:e|E)igenvalues of the unoccupied subspace spin(.*?\d{8}[\.\-\d\s]*)',
            rundata, re.DOTALL)
        levelregex = '\-?\d\.\d{8}'

        fermidatas = re.findall('Fermi Energy.*', rundata)
        fermiregex = '\-?\d\.\d{6}'

        self.energylevels = []
        self.spins = []
        for i in range(0, len(occupieddatas)):
            # If we have unoccupied levels...
            if len(unoccupieddatas) == len(occupieddatas):
                E = np.array(
                    re.findall(
                        levelregex, occupieddatas[i] + unoccupieddatas[i]),
                    dtype=float)
            else:
                E = np.array(re.findall(
                    levelregex, occupieddatas[i]), dtype=float)
            E = E * atc.Ha / atc.eV
            fermi = float(re.search(fermiregex, fermidatas[i]).group(0))

            spin = str(re.search('\s*(\d)', occupieddatas[i]).group(1))
            spin = int(spin) - 1

            levels = fu.EnergyLevels(energies=E, fermi=fermi)
            self.energylevels.append(levels)
            self.spins.append(spin)

    def read_from_pdos(self, fname):
        """Reads Spectrum from projected density of states"""
        # format: id energy occupation s py pz px d-2 d-1 d0 d+1 d+2 f-3 f-2 f-1 f0 f+1 f+2 f+3
        A = np.genfromtxt(fname, skip_header=2)

        energies = A[:, 1] * atc.Ha / atc.eV
        occupations = np.array(A[:, 2], dtype=int)
        weights = np.sum(A[:, 3:], axis=1)

        # needed only for fermi energy
        s = open(fname, 'r').readline()
        print(s)
        fermiregex = 'E\(Fermi\) =\s+([-\d\.]+) a.u.'
        match = re.search(fermiregex, s)
        if not match:
            raise ValueError(
                "Unable to parse Fermi energy from pdos file. Setting to zero.")
            fermi = 0
        else:
            fermi = float(match.group(1)) * atc.Ha / atc.eV

        self.energylevels = [fu.EnergyLevels(energies=energies,
                                             occupations=occupations, weights=weights, fermi=fermi)]
        self.spins = [0]


class WfnCube(cube.Cube):
    """Gaussian cube file written by CP2K

    CP2K writes the index of level and spin into the
    comment line of the cube file
    """

    def __init__(self, title=None, comment=None, origin=None, atoms=None,
                 data=None, spin=None, wfn=None, energy=None, occupation=None):
        """Standard constructor, all parameters default to None.

        energy and occupation are not stored in the cube file,
        but can be assigned by linking the cube file with the 
        output from the calculation.
        """
        super(WfnCube, self).__init__(title, comment, origin, atoms, data)
        self.spin = spin
        self.wfn = wfn
        self.energy = energy
        self.occupation = occupation

    @classmethod
    def from_file(cls, fname, read_data=False):
        """Creates Cube from cube file"""
        tmp = WfnCube()
        tmp.read_cube_file(fname, read_data=read_data)
        return tmp

    def read_cube_file(self, fname, read_data=False, v=1):
        """Reads header and/or data of cube file"""
        super(WfnCube, self).read_cube_file(fname, read_data, v)

        # CP2K stores information on the level/spin index
        # in the comment line
        commentregex = 'WAVEFUNCTION\s+(\d+)\s+spin\s+(\d+)'
        match = re.search(commentregex, self.comment)
        self.wfn = int(match.group(1))
        self.spin = int(match.group(2))
