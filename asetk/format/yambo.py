"""Classes for use with Yambo

Representation of a spectrum.
Main functionality is to read from Yambo output, o.qp files
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

    @property
    def nk(self):
        return len(self.energylevels)


class Spectrum(object):

    """A Spectrum holds the data belonging to all spins"""

    def __init__(self, energylevels=None):
        """Set up spectrum from a list of EnergyLevels."""
        self.dispersions = [ Dispersion(energylevels) ]

    @classmethod
    def from_output(cls, fname, mode='QP'):
        """Creates Spectrum from Yambo output file"""
        tmp = Spectrum()
        tmp.read_from_output(fname, mode)
        return tmp

    @classmethod
    def from_qp(cls, fname=None, mode='QP'):
        """Creates Spectrum from Yambo o.qp file"""
        tmp = Spectrum()
        tmp.read_from_qp(fname, mode)
        return tmp

    @classmethod
    def from_netcdf_db(cls, fname=None, mode='QP'):
        """Creates Spectrum from Yambo netcdf database"""
        tmp = Spectrum()
        tmp.read_from_netcdf_db(fname, mode=mode)
        return tmp

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

    @property
    def nspin(self):
        return len(self.dispersions)

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

    def read_from_output(self, fname, mode=None):
        s = open(fname, 'r').read()

        floatregex = '-?\d+\.\d+'
        lineregex='[^\r\n]*\r?\n'
        #blanklineregex='(?:^\s*$)'

        if mode == 'DFT' or mode == None:
            kptregex = 'X\* K.*?: ({f})\s*({f})\s*({f}).*?weight\s*({f}){l}(.*?)[\*\[]'\
                        .format(f=floatregex,l=lineregex)
            fermiregex='Fermi Level.*?:(\s*[\-\d\.]+)'
        elif mode == 'QP':
            kptregex = 'Q?P \[eV\].*?:\s*({f})\s+({f})\s+({f})(.*?)[Q\[]'\
                        .format(f=floatregex)


        self.spins=[]
        self.dispersions=[]

        # No spin for the moment, but shouldn't be too difficult to extend
        for spin in [0]:

            disp = Dispersion()
            matches=re.findall(kptregex, s, re.DOTALL)

            if mode == 'DFT' or mode == None:

                fermi = float(re.search(fermiregex, s).group(1))
                energylevels = []
                kvectors = []
                weights = []

                for match in matches:
                    kx, ky, kz, weight, ldata  = match
                    kvectors.append( np.array([kx, ky, kz], dtype=float) )
                    weights.append( float(weight) )
                    energies = re.findall('({f})'.format(f=floatregex), ldata, re.DOTALL)


                    energies = np.array(energies, dtype=float)
                    levels = fu.EnergyLevels(energies=energies,occupations=None, fermi=fermi)
                    energylevels.append(levels)
                 
                disp = Dispersion(energylevels=energylevels, kvectors=kvectors, weights=weights)

            elif mode == 'QP':

                energylevels = []
                kvectors = []

                for match in matches:
                    kx, ky, kz, ldata  = match
                    kvectors.append( np.array([kx, ky, kz], dtype=float) )
                    energies = re.findall('E=\s*({f})'.format(f=floatregex), ldata, re.DOTALL)

                    energies = np.array(energies, dtype=float)
                    levels = fu.EnergyLevels(energies=energies)
                    energylevels.append(levels)
                 
                disp = Dispersion(energylevels=energylevels, kvectors=kvectors)

            self.dispersions.append(disp)
            self.spins.append(spin)

        
    def read_from_qp(self, fname="o.qp", ihomo=None):
        """Read from o.qp output (has more digits than in report.

           Anyhow, the proper way would be to read the database"""
        s = open(fname, 'r').read()

        data = np.genfromtxt(fname, dtype=float)

        energies = data[:,2] + data[:,3]

        # setting HOMO to zero
        if ihomo:
            energies -= energies[ihomo]

             
        self.spins=[]
        self.dispersions=[]

        # No spin for the moment, but shouldn't be too difficult to extend
        for spin in [0]:
            levels = fu.EnergyLevels(energies=energies,occupations=None)
            disp = Dispersion(energylevels=[levels], kvectors = [ (0,0,0) ] )


        self.dispersions.append(disp)
        self.spins.append(spin)

    def read_from_netcdf_db(self, fname="ndb.QP", mode="QP"):
        """Read from netCDF database
        
        requires netCDF4 python module"""

        from netCDF4 import Dataset
        f = Dataset(fname, 'r')
        SPIN_VARS = f.variables['SPIN_VARS'][:]
        QP_kpts =  f.variables['QP_kpts'][:]
        QP_table = f.variables['QP_table'][:]
        QP_E_Eo_Z =  f.variables['QP_E_Eo_Z'][:]
        f.close()
        
        nspin = len(SPIN_VARS)

        nk = QP_kpts.shape[1]
        kpts = [ QP_kpts[:,ik] for ik in range(nk) ]

        ibnds, dum, iks, ispins = QP_table
        nbnd = len(ibnds) / (nspin * nk)

        if mode == "QP":
            iener = 0
        elif mode == "DFT":
            iener = 1
        else:
            print("Error: Did not recognize mode '{}'.".format(mode))

        self.spins=[]
        self.dispersions=[]
        for ispin in range(nspin):
            is_spin = np.where(ispins == SPIN_VARS[ispin])[0]

            energylevels = []
            kvectors = []
            for ik in range(nk):
                k = kpts[ik]

                is_k = np.where(iks == ik+1)[0]
                # still need to figure out the first index
                # is it real vs. complex?
                e = QP_E_Eo_Z[0, np.intersect1d(is_spin,is_k), iener] * atc.Ha / atc.eV
                levels = fu.EnergyLevels(energies=e,occupations=None)

                kvectors.append(k)
                energylevels.append(levels)

            disp = Dispersion(energylevels=energylevels, kvectors = kvectors)

            self.dispersions.append(disp)
            self.spins.append(ispin)

        ## setting HOMO to zero
        #if ihomo:
        #    energies -= energies[ihomo]

