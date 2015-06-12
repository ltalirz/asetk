"""Classes for use with BerkeleyGW

Representation of a spectrum.
Main functionality is to read from BerkeleyGW output, eqp.dat files
and also hdf5 databases.
"""

import re
import copy  as cp
import numpy as np
import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as atc
from . import cube


class Spectrum(object):

    """A Spectrum holds the data belonging to all spins"""

    def __init__(self, dispersions=None, spins=None):
        """Set up spectrum from list of dispersions."""
        self.dispersions = dispersions
        self.spins = spins

    @classmethod
    def from_eqp(cls, fname=None, mode=None):
        """Creates Spectrum from eqp.dat file"""
        tmp = Spectrum()
        tmp.read_from_eqp(fname, mode)
        return tmp

    @classmethod
    def from_log(cls, fname=None, mode=None):
        """Create Spectrum from .log file
        
        These files are written by sigma.x"""
        tmp = Spectrum()
        tmp.read_from_log(fname, mode=mode)
        return tmp

    @classmethod
    def from_hdf5_db(cls, fname=None, mode=None):
        """Create Spectrum from eps0mat.h5 file
        
        These files are written by epsilon.x"""
        tmp = Spectrum()
        tmp.read_from_hdf5_db(fname, mode=mode)
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
    def fermi(self):
        """Returns Fermi energy."""
        fermis = [d.fermi for d in self.dispersions]

        fermi = np.unique(fermis)

        if len(fermi) == 1:
            return fermi[0]
        elif len(fermi) != 1:
            print("There are Fermi energies {}".format(fermis))
            print("Using the mean {}".format(np.mean(fermis)))
            return np.mean(fermis)


    @property
    def occupations(self):
        """Returns list of level occupations of all spins."""
        os = []
        for disp in self.dispersions:
            os = os + disp.occupations
        return os

    @property
    def nbnd(self):
        nbnds = [d.nbnd for d in self.dispersions]
        nbnd = np.unique(nbnds)

        if len( np.unique(nbnd) ) != 1:
            print("Warning: spins have different numer of bands {}"\
                   .format(nbnd))
        return nbnd[0]

    @property
    def nkpt(self):
        nkpts = [d.nkpt for d in self.dispersions]
        nkpt = np.unique(nkpts)

        if len( np.unique(nkpt) ) != 1:
            print("Warning: spins have different numer of k-points {}"\
                   .format(nkpt))
        return nkpt[0]

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
        return self.dispersions[index]

    def read_from_eqp(self, fname="sigma_hp.log", ihomo=None):
        """Read from eqp.dat files generated with eqp.py script
        """
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
                     
            disp = fu.Dispersion(energylevels=energylevels, kvectors=kvectors)
            self.dispersions.append(disp)
            self.spins.append(spin)

    def read_from_log(self, fname="sigma_hp.log", mode=None, ihomo=None):
        """Read from sigma_hp.log files generated by sigma.x
        """
        s = open(fname, 'r')
        string = s.read()
        s.close()

        # get number of bands
        match = re.search('band_index\s+(\d+)\s+(\d+)', string)
        ibnd_min = int(match.group(1))
        ibnd_max = int(match.group(2))
        nbnd = ibnd_max - ibnd_min + 1

        # get number of spins
        smatches = re.findall('spin \= (\d+)', string)
        spins = np.array(smatches, dtype=int)
        self.spins = np.unique(spins) - 1 # want to start from 0
        nspin = len(self.spins)

        kregex = '  k \=.*?(?=\n{3}|   k \=)'
        kmatches = re.findall(kregex, string, re.DOTALL)

        if mode == 'QP' or mode == None:
            if re.search("ch'",kmatches[0]) is not None:
                calc_mode = 'static_remainder'
            else:
                calc_mode = 'partial_sum'
        else:
            raise ValueError("Unknown mode {}".format(mode))

        spins = [[] for s in self.spins]
        for kmatch in kmatches:
            lines = kmatch.splitlines()

            dum, dum, kx, ky, kz, dum, dum, ik, dum, dum, ispin = lines[0].split()
            ispin = int(ispin)-1
            kvector = np.array([kx, ky, kz], dtype=float)

            energies = []
            for line in lines[3:3+nbnd]:
                if calc_mode == 'static_remainder':
                    n, elda, ecor, x, sx_x, ch ,sig, vxc, eqp0, eqp1, \
                        ch_p, sig_p, eqp0_p, eqp1_p, Znk = line.split()
                else:
                    n, elda, ecor, x, sx_x, ch ,sig, vxc, eqp0, eqp1,\
                        Znk = line.split()
                energies.append(float(eqp1))

            levels = fu.EnergyLevels(energies=np.array(energies, dtype=float))

            kpt = fu.KPoint(kvector=kvector, energylevels=levels)
            spins[ispin].append(kpt)
                     
        self.dispersions=[]
        for s in spins:
            disp = fu.Dispersion(kpoints=s)
            self.dispersions.append(disp)


    def read_from_hdf5_db(self, fname="eps0mat.h5", mode="QP"):
        """Read from HDF5 database
        
        requires netCDF4 python module"""

        from netCDF4 import Dataset

        f = Dataset(fname, 'r')
        #print(f.groups)
        kpoints_grp = f.groups['mf_header'].groups['kpoints']
        #kpoints_grid = kpoints_grp.variables['kgrid'][:]
        kpoints_rk = kpoints_grp.variables['rk'][:]
        #kpoints_nrk = kpoints_grp.variables['nrk'][:]
        kpoints_el = kpoints_grp.variables['el'][:]
        kpoints_nspin = kpoints_grp.variables['nspin'][:]
        kpoints_mnband = kpoints_grp.variables['mnband'][:]
        f.close()

        nk = kpoints_rk.shape[0]
        nbnd = kpoints_mnband
        nspin = kpoints_nspin

        #if mode == "QP":
        #    iener = 0
        #elif mode == "DFT":
        #    iener = 1
        #else:
        #    print("Error: Did not recognize mode '{}'.".format(mode))

        self.spins=[]
        self.dispersions=[]

        for ispin in range(nspin):

            self.spins.append(ispin)

            kpoints = []
            for ik in range(nk):
                k = kpoints_rk[ik]
                e = kpoints_el[ispin, ik, :] * atc.Ry / atc.eV
                levels = fu.EnergyLevels(energies=e,occupations=None)

                kpt = fu.KPoint(kvector=k, energylevels=levels)
                kpoints.append(kpt)

            disp = fu.Dispersion(kpoints=kpoints)
            self.dispersions.append(disp)

        ## setting HOMO to zero
        #if ihomo:
        #    energies -= energies[ihomo]

