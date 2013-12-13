"""Shorthand notations for physical constants

This module is using scipy.constants.physical_constants.
All constants are in SI units in order to allow for arbitrary
conversions between them.
"""

from scipy.constants import physical_constants

# Energies
eV = physical_constants['electron volt-joule relationship'][0]
Ha = physical_constants['Hartree energy'][0]
Ry = Ha / 2.0

# Distances
a0 = physical_constants['Bohr radius'][0]
Angstrom = 1e-10

