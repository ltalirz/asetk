"""Classes for use with Quantum ESPRESSO

Representation of a spectrum
"""

import re
import os
import copy  as cp
import numpy as np

import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as atc

#class Dispersion(object):
#    """Holds k-points belonging to one spin"""
#    def __init__(self, kpoints=None, kvectors=None):
#        if kpoints is None:
#            self.kpoints = []
#            self.kvectors = []
#        else:
#            self.kpoints = kpoints
#            self.kvectors = kvectors
#
#    def addkpoint(self, kpoint, kvector):
#        self.kpoints.append(kpoint)
#        self.kvectors.append(kvector)
#
#    def energylevels(self):
#        s = fu.EnergyLevels()
#        for l in self.kpoints:
#            s += l
#        return s
#
#    def mergekpoints(self):
#        self.kpoints = energylevels(self.kpoints)
#
#    def __iadd__(self, s):
#        """Merging two spins"""
#        if len(self.kpoints) != len(s.kpoints):
#            print("Unable to merge due to different number of kpoints")
#        for i in range(len(self.kpoints)):
#            self.kpoints[i] += s.kpoints[i]
#        return self
#
#    def shift(self, e):
#        for l in self.kpoints:
#            l.shift(e)
#
#    @property
#    def nbnd(self):
#        nbnds = [len(k.levels) for k in self.kpoints]
#        nbnd = np.unique(nbnds)
#
#        if len( np.unique(nbnd) ) != 1:
#            print("Warning: k-points have different numer of bands {}"\
#                   .format(nbnd))
#        return nbnd[0]
#
#    @property
#    def nk(self):
#        return len(self.kpoints)



class Spectrum(object):
    """A collection of dispersions, grouped by spin"""

    def __init__(self, dispersions=None, spins=None):
        """Set up spectrum from list of dispersions."""
        self.dispersions = dispersions
        self.spins = spins

    @classmethod
    def from_save(cls, fname):
        """Creates Spectrum from QE save directory"""
        tmp = Spectrum()
        tmp.read_from_save(fname)
        return tmp

    @classmethod
    def from_output(cls, fname):
        """Creates Spectrum from QE output"""
        tmp = Spectrum()
        tmp.read_from_output(fname)
        return tmp

    @property
    def energies(self):
        """Returns list of all energy levels of all spins."""
        list = [d.energies for d in self.dispersions]
        return np.concatenate(list)

    @property
    def occupations(self):
        """Returns list of level occupations of all spins."""
        os = np.array([])
        for d in self.dispersions:
            os = np.concatenate( (os, d.occupations))
        return os

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
        self.dispersions = [ d.copy() for d in spectrum.dispersions ]
        self.spins = cp.copy(spectrum.spins)

    def shift(self, de):
        for d in self.dispersions:
            d.shift(de)

    def __str__(self):
        text  = "Spectrum containing {} spins\n".format(len(self.dispersions))
        for i in range(len(self.dispersions)):
            d = self.dispersions[i]
            s = self.spins[i]
            text += 'spin {} : {}\n'.format(s+1, d.__str__())
        return text

    def __getitem__(self, index):
        return self.dispersions[index]

    def read_from_save(self, prefix):
        """Reads Spectrum from QE save directory"""

        savedir = prefix + '.save'
        if not os.path.exists(savedir):
            raise IOError("Directory {s} not found.".format(s=savedir))
        os.chdir(savedir)

        dataxml = open('data-file.xml', 'r').read()

        nspinregex = '<NUMBER_OF_SPIN_COMPONENTS.*?>\s*(\d+)'
        nspin = int( re.search(nspinregex, dataxml, re.DOTALL).group(1) )

        # should be able to match scientific and non-scientific notation
        floatregex = '-?\d+\.\d+(?:[Ee][+\-]?\d+)?'

        # get fermi energy
        fermiregex = '<FERMI_ENERGY.*?>\s*({f})'.format(f=floatregex)
        fermi = float( re.search(fermiregex, dataxml, re.DOTALL).group(1) )
        fermi *= atc.Ha / atc.eV

        #get lattice parameter
        alatregex = '<LATTICE_PARAMETER.*?>\s*({f})'.format(f=floatregex)
        alat = float( re.search(alatregex, dataxml, re.DOTALL).group(1) )
        alat *= atc.a0  


        # get k-points
        kptregex = '<K-POINT\.(\d+)\s+XYZ=\"({f}) ({f}) ({f})\"\s+WEIGHT=\"({f})\"/>'\
                .format(f=floatregex)
        kptdatas = re.findall(kptregex, dataxml)

        self.dispersions = []
        self.spins = []
        for spin in range(nspin):
            dispersion = fu.Dispersion()

            for kpt in kptdatas:
                kindex = int(kpt[0])
                kvec   = np.array([ kpt[1], kpt[2], kpt[3] ], dtype=float)
                #kvec  *= 2*np.pi / alat

                kdir = 'K{k:05d}'.format(k=kindex)
                if not os.path.exists(kdir):
                    raise IOError("Directory {s} not found.".format(s=kdir))

                # Read energy levels
                os.chdir(kdir)

                # get correct file name
                if nspin == 1:
                    eigf = 'eigenval.xml'
                elif nspin == 2:
                    eigf = 'eigenval{}.xml'.format(spin+1)
                else:
                    print("Error: Can only handle nspin=1, 2")

                if not os.path.exists(eigf):
                    print("Error: Cannot find file {}".format(eigf))

                eigenvalxml = open(eigf, 'r').read()

                eigregex  = '<EIGENVALUES.*?>(.*?)<'
                eigstring = re.search(eigregex, eigenvalxml, re.DOTALL).group(1)

                levelregex  = '\-?\d+.*'
                levelregex  = floatregex
                levels = np.array(re.findall(levelregex, eigstring), dtype = float)
                levels *= atc.Ha / atc.eV

                kpt = fu.KPoint(
                        energylevels=fu.EnergyLevels(energies=levels, fermi=fermi),
                        kvector=kvec
                        )
                dispersion.kpoints.append(kpt)

                os.chdir('..')

            self.dispersions.append(dispersion)    
            self.spins.append(spin)    

        os.chdir('..')

    def read_from_output(self, prefix):
        return 0;


