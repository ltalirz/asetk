"""Classes for use with BerkeleyGW

Representation of a spectrum.
Main functionality is to read from BerkeleyGW output, eqp.dat files
and also netcdf databases.
"""

import re
import copy  as cp
import numpy as np
import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as atc
from . import cube

class Dispersion:
    """A Dispersion holds the k-points belonging to one spin"""

    def __init__(self, energylevels=None, kvectors=None, weights=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.__energylevels = energylevels
        self.__kvectors = kvectors
        self.__weights = weights

    @property
    def energylevels(self):
        """Returns energylevelsi of all k-points."""
        return self.__energylevels

    @property
    def kvectors(self):
        return self.__kvectors

    @property
    def weights(self):
        return self.__weights

    @property
    def energies(self):
        """Returns list of energy levels of all k-points."""
        list = [el.energies for el in self.__energylevels]
        return np.concatenate(list)

    @property
    def occupations(self):
        """Returns list of level occupations of all k-points."""
        os = []
        for el in self.__energylevels:
            os = os + list(el.occupations)
        return os

    def copy(self, dispersion):
        """Performs deep copy of dispersion."""
        self.__energylevels = [ el.copy() for el in dispersion.__energylevels ]
        self.__kvectors = cp.copy(spectrum.__kvectors)
        self.__weights = cp.copy(spectrum.__weights)

    def shift(self, de):
        for levels in self.__energylevels:
            levels.shift(de)

    def __str__(self):
        text  = "Dispersion containing {} k-points\n".format(len(self.__energylevels))
        for i in range(len(self.__energylevels)):
            e = self.__energylevels[i]
            k = self.__kvectors[i]
            text += 'k = ({:6.3f}, {:6.3f}, {:6.3f})'.format(k[0], k[1], k[2])
            if self.__weights:
                w = self.__weights[i]
                text += ', w = {}'.format(w)
            text += ' : {}\n'.format(e.__str__())
        return text

    def __getitem__(self, index):
        return self.__energylevels[index]


class Spectrum(object):

    """A Spectrum holds the data belonging to all spins"""

    def __init__(self, energylevels=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.dispersions = [ Dispersion(energylevels) ]

    @classmethod
    def from_output(cls, fname, mode=None):
        """Creates Spectrum from Yambo output file"""
        tmp = Spectrum()
        tmp.read_from_output(fname, mode)
        return tmp

    @classmethod
    def from_qp(cls, fname=None, mode=None):
        """Creates Spectrum from Yambo o.qp file"""
        tmp = Spectrum()
        tmp.read_from_qp(fname, mode)
        return tmp

##@classmethod
#    def from_netcdf_db(cls, fname=None, mode=None):
#        """Creates Spectrum from Yambo netcdf database"""
#        tmp = Spectrum()
#        tmp.read_from_netcdf_db(fname, mode=mode)
#        return tmp

    @property
    def energies(self):
        """Returns list of energies e[ispin][ibnd]."""
        list = [disp.energies for disp in self.dispersions]
        return list

    @property
    def energylevels(self):
        """Returns list of Energylevels l[ispin][ibnd]."""
        list = []
        for d in self.dispersions:
            sum = fu.Energylevels()
            for el in d.energylevels:
                sum += el
            list.append(sum)
        return list

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

    def read_from_qp(self, fname="o.qp", ihomo=None):
        """Read from o.qp output (has more digits than in report.

           Anyhow, the proper way would be to read the database"""
        s = open(fname, 'r')

        self.spins=[]
        self.dispersions=[]

        # No spin for the moment, but shouldn't be too difficult to extend
        for spin in [0]:

            energylevels = []
            kvectors = []

            while True:
                # 1st line of each block contains k-vector and number of bands
                line = s.readline()
                if not line:
                    break
                kx, ky, kz, nbnd = line.split()
                #print(kx, ky, kz, nbnd)
                kvectors.append( np.array([kx, ky, kz], dtype=float) )

                energies = []
                for i in range(int(nbnd)):
                    ispin, ibnd, edft, eqp = s.readline().split()
                    energies.append(eqp)

                levels = fu.EnergyLevels(energies=np.array(energies, dtype=float))
                energylevels.append(levels)
                     
            disp = Dispersion(energylevels=energylevels, kvectors=kvectors)
            self.dispersions.append(disp)
            self.spins.append(spin)


    def read_from_netcdf_db(self, fname="ndb.QP", mode="QP"):
        """Read from netCDF database
        
        requires netCDF4 python module"""

        from netCDF4 import Dataset
        f = Dataset(fname, 'r')

        QP_kpts =  f.variables['QP_kpts'][:]
        QP_table = f.variables['QP_table'][:]
        QP_E_Eo_Z =  f.variables['QP_E_Eo_Z'][:]

        nk = QP_kpts.shape[1]
        nbnd = QP_table.shape[1] / nk
        nspin = QP_E_Eo_Z.shape[0]
        if mode == "QP":
            iener = 0
        elif mode == "DFT":
            iener = 1
        else:
            print("Error: Did not recognize mode '{}'.".format(mode))

        self.spins=[]
        self.dispersions=[]

        for ispin in range(nspin):

            energylevels = []
            kvectors = []
            for ik in range(nk):
                k = QP_kpts[:, ik]
                e = QP_E_Eo_Z[0, ik*nbnd : (ik+1)*nbnd, iener] * atc.Ha / atc.eV
                levels = fu.EnergyLevels(energies=e,occupations=None)

                energylevels.append(levels)
                kvectors.append(k)

            disp = Dispersion(energylevels=energylevels, kvectors = kvectors)

            self.dispersions.append(disp)
            self.spins.append(ispin)

        ## setting HOMO to zero
        #if ihomo:
        #    energies -= energies[ihomo]

