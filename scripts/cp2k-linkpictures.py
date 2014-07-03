#!/usr/bin/env python
# Link png pictures of orbitals to the corresponding energy, writes file 'list.dat'
# Usage: cp2k-linkpictures.py cp2k.out "*.png"
import re
import argparse
import glob
import os

import asetk.format.cp2k as cp2k


# Define command line parser
parser = argparse.ArgumentParser(
    description='Outputs a list linking files to their respective energy level')
parser.add_argument('--version', action='version', version='%(prog)s 22.01.2014')
parser.add_argument(
    'out',
    metavar='WILDCARD', 
    help='CP2K output containing energy levels')
parser.add_argument(
    'pics',
    metavar='WILDCARD', 
    help='Some files belonging to an energy level, format: PROJ-WFN_00329_2-1_0.png')

args = parser.parse_args()
expanded = glob.glob(args.pics)

spectrum = cp2k.Spectrum.from_mo(args.out)

outputfile = 'list.dat'
f=open(outputfile,'w')

for filename in expanded:
    spin =  int(re.search('_(\d)\-', filename).group(1))
    level =  int(re.search('WFN_(\d+)', filename).group(1))

    print >>f, "{spin}\t{l}\t{e}\t{f}".format(spin=spin, l=level, e=spectrum.energylevels[spin-1].energies[level-1], f=filename)
    #print >>f, "%(spin)s\t%(nLevel)s\t%(level)s\t%(file)s" % dict(spin=spin, nLevel=nLevel, level=energies[spin-1][nLevel-1], file=filename)
    
f.close()
