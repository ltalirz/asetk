"""Classes for use with Quantum ESPRESSO

Representation of a spectrum
"""

import re
import os
import copy  as cp
import numpy as np
from string import digits

import asetk.atomistic.fundamental as fu
import asetk.atomistic.constants as atc
from . import cube

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



class QECube:
    """Intermediate cube file format written by pp.x

    These files contain (squared) Kohn-Sham wave functions (?)
    in a text-based format very similar to the Gaussian cube format.
    
    Format specification
    
    LINE   FORMAT      CONTENTS
    ===============================================================
     1     A           TITLE
     2     8I8         NX NY NZ NX NY NZ #ATOMS #SPECIES
     3     I8,6F16.8   IBRAV CELLDM(1:6)
     4-6   3F16.8      CELL VECTORS IN ALAT UNITS, ONLY IF IBRAV=0
                       (OTHERWISE CONTINUE WITH LINE 7)
     7     3F16.8,I8   ???
     #SPECIES LINES OF ATOMIC SPECIES INFO:
     ...   I4,S4,6.2F  INDEX, LABEL, #VALENCE ELECTRONS OF PSEUDO
     #ATOMS LINES OF ATOM COORDINATES:
     ...   I5,3F12.6,I4 ATOM INDEX, X, Y, Z [ALAT UNITS], SPECIES INDEX
     REST: 5E17.9      CUBE DATA (WITH X INCREMENT MOVING FASTEST, THEN
                       Y AND THEN Z)
    
    ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
    """

    def __init__(self, filename=None, title=None, cell=None, atoms=None,
            data=None):
        """Standard constructor, all parameters default to None"""

        self.filename = filename
        self.title = title
        self.cell = cell
        self.atoms = atoms
        self.data = data
        self._shape = None  # stores shape, if grid data isn't read

    @classmethod
    def from_file(cls, fname, read_data=False):
        """Creates Cube from cube file"""
        tmp = QECube()
        tmp.read_qe_cube_file(fname, read_data=read_data)
        return tmp

    def read_qe_cube_file(self, fname, read_data=False, v=1):
        """Reads header and/or data of cube file"""
        #super(WfnCube, self).read_cube_file(fname, read_data, v)
        self.filename = fname
        b2A = atc.a0 / atc.Angstrom

        f = open(fname, 'r')
        readline = f.readline

        # line 1
        self.title = readline()
#       self.comment = readline()

#       axes = [0, 1, 2]
        # line 2
        line = readline().split()
        nx, ny, nz, nxs, nys, nzs, natoms, nspecies = np.array(line, dtype=int)
        shape = np.array([nx, ny, nz], dtype=int)
        self._shape = shape

        # line 3
        line = readline().split()
        ibrav = int(line[0])
        alat = float(line[1])
        celldm = np.array(line[2:],dtype=float)
        if ibrav > 0:
            raise ValueError("ibrav > 0 not yet implemented.")

        # lines 4-6
        cell = np.empty((3,3))
        for i in range(3):
            x, y, z = [float(s) for s in readline().split()]
            cell[i] = np.array([x,y,z], dtype=float)
        cell *= b2A * alat

        # line 7
        line = readline().split()

        # species
        species = np.empty(nspecies, dtype=str)
        for i in range(nspecies):
            sindex, symbol, valence_electrons = readline().split()
            # removing any digits that may be part of the symbol
            # such as C1, C2, ...
            species[i] = symbol.translate(None, digits)

        # atoms
        at_positions = np.empty((natoms, 3))
        at_symbols = np.empty(natoms, dtype=str)
        for i in range(natoms):
            line = readline().split()
            at_positions[i] = [float(s) for s in line[1:4]]   
            at_symbols[i] = species[int(line[4])-1]

        at_positions *= b2A * alat
        pbc = [True, True, True]
        self.atoms = fu.Atoms(symbols=at_symbols, positions=at_positions, 
                cell=cell, pbc=pbc)

        if read_data:
            # Note:
            # This is already ~1.7x faster than ASE's version.
            # However, parsing still dominates for reasonable disk speeds
            # (parsing time = 8x reading time on 480 MB/s SSD)
            # In order to parse quickly, use read_csv from the pandas module.

            # read() pretty much maxes out the disk read speed.
            # split() takes a considerable amount of time.
            # The conversion to float is even more expensive.
            self.data = np.array(f.read().split(), dtype=float)

            # In QE's format, the fastest index is x, then y, then z
            self.data = self.data.reshape(shape[::-1])
            self.data = self.data.swapaxes(0,2)

        f.close()

    def write_cube_file(self, fname=None):
        """Writes object to Gaussian cube file

        """
	tmp = self.to_cube()
        tmp.write_cube_file()

    def to_cube(self):
        """Converts object to Gaussian cube object

        """
        tmp = cube.Cube(filename=fname, title=self.title, 
                comment="Converted from QE intermediate cube format\n",
                origin = np.array([0,0,0]),
                atoms = self.atoms,
                data = self.data)

	return tmp
