#!/usr/bin/env python
import numpy as np
import scipy.special as scsp
import argparse
import asetk.format.cp2k as cp2k
import asetk.util.progressbar as progressbar
import asetk.atomistic.constants as constants
import os.path

# Define command line parser
parser = argparse.ArgumentParser(
    description='Performs Scanning Tunneling Spectroscopy Simulation.')
parser.add_argument('--version', action='version', version='%(prog)s 13.12.2013')
parser.add_argument(
    'cubes',
    nargs='+',
    metavar='WILDCARD',
    help='Gaussian cube files containing the Kohn-Sham orbitals psi_i \
          (or |psi_i|^2, see --psisquared).')
parser.add_argument(
    '--psisquared',
    dest='psisquared',
    action='store_true',
    default=False,
    help='If specified, the cubes are assumed to contain |psi_i|^2.')
parser.add_argument(
    '--levelsfile',
    metavar='FILENAME',
    required=True,
    help='File containing the energy levels. Can be either CP2K output\
          or .MOLog file.')
parser.add_argument(
    '--outfile',
    metavar='FILENAME',
    default='sts.cube',
    help='Name of cube file for STS simulation.')
parser.add_argument(
    '--emin',
    metavar='ENERGY',
    default=-3.0,
    type=float,
    help='Minimum bias for STS [eV].')
parser.add_argument(
    '--emax',
    metavar='ENERGY',
    default=+3.0,
    type=float,
    help='Maximum bias for STS [eV].')
parser.add_argument(
    '--de',
    metavar='ENERGY',
    default=+0.01,
    type=float,
    help='Bias step for STS [eV].')
parser.add_argument(
    '--eref',
    metavar='ENERGY',
    default=None,
    type=float,
    help='By default, Fermi is taken as zero-energy reference. With this option\
          it is possible to define a different reference')
parser.add_argument(
    '--FWHM',
    metavar='ENERGY',
    default=+0.100,
    type=float,
    help='Full-width half-maximum of broadening.\n\
          FWHM = 2.355 \sigma in Gaussian broadening. \
          FWHM = \Gamma in Lorentzian broadening.')
parser.add_argument(
    '--bmethod',
    metavar='KEYWORD',
    default='Gaussian',
    help='The broadening method applied.\n\
          Can be \'Gaussian\' or \'Lorentzian\'.')
parser.add_argument(
    '--bepsilon',
    metavar='FLOAT',
    default=1e-5,
    type=float,
    help='Convolution is performed with broadening function in range \n\
          [E-Eb, E+Eb] where Eb is determined such that the integrated \n\
          weight of the broadening function outside of [-Eb, Eb] is \n\
          below bepsilon.')
parser.add_argument(
    '--height',
    metavar='DISTANCE',
    default=2.5,
    type=float,
    help='The height [Angstroms] above the topmost atom, where to extract the plane.')
parser.add_argument(
    '--normalize',
    dest='normalize',
    action='store_true',
    default=False,
    help='If provided, the STS intensity is normalized to 1.')
parser.add_argument(
    '--print-weights',
    metavar='FILENAME',
    dest='print_weights',
    default=None,
    help='If specified, the (relative) weights   \
            w_i = \int |\psi_i^2(x,y,z_0)| dx dy \
          of all states are printed to file.')

args = parser.parse_args()

# Prepare broadening functions for later convolution
# quantile function quantile(y) = x is defined such that
# the integral of the probability density from -\infty to x equals y.
if args.bmethod == 'Gaussian':
    sigma = args.FWHM / np.log(8 * np.sqrt(2))
    quantile = lambda y: np.sqrt(2) * scsp.erfinv(2.0*y -1.0) * sigma
    broadening = lambda x: 1/(sigma * np.sqrt(2*np.pi)) \
                         * np.exp( - x**2 / (2 * sigma**2) )
elif args.bmethod == 'Lorentzian':
    gamma = args.FWHM * 0.5
    quantile = lambda y: np.tan(np.pi * (y - 0.5)) * gamma
    broadening = lambda x: 1/np.pi * gamma / (x**2 + gamma**2)
else:
    print("Error: No broadening method recognized.")

eb = -quantile(args.bepsilon * 0.5)
print("Using window [x{:.3f}eV, x+{:.3f}eV] for level broadening."\
      .format(-eb, eb))
print("This window contains {:.5f} % of the total weight of {}.\n"\
      .format((1-args.bepsilon)*100, args.bmethod))

lfname, lfext = os.path.splitext(args.levelsfile)
if lfext == '.MOLog':
    spectrum = cp2k.Spectrum.from_mo(args.levelsfile)
else:
    spectrum = cp2k.Spectrum.from_output(args.levelsfile)


print("Read spectrum from {f}".format(f=args.levelsfile))
print(spectrum)

if args.eref is not None:
    spectrum.shift(-args.eref)
    print("Taking {s} eV as zero energy reference.".format(s=args.eref))
else:
    for spin, levels in zip(spectrum.spins, spectrum.energylevels):
        levels.shift(-levels.fermi)
    print("Fermi energy is taken as zero energy reference.")

# Reading headers of cube files
cubes = []
print("\nReading headers of {n} cube files".format(n=len(args.cubes)))
bar = progressbar.ProgressBar(niter=len(args.cubes))

for fname in args.cubes:
    cubes.append( cp2k.WfnCube.from_file(fname) )
    bar.iterate()
print("")

# Connecting energies with required cube files
required_cubes = []
for spin, levels in zip(spectrum.spins, spectrum.energylevels):
    for index, l in enumerate(levels.levels):
        e = l.energy
        o = l.occupation
        # If we need the cube file for this level..
        if e >= args.emin - eb and \
           e <= args.emax + eb:

               found = False
               for cube in cubes:
                   #print "ind {}  {} spin {} {}".format(cube.wfn,index,cube.spin,spin)
                   if cube.wfn == index+1 and cube.spin == spin + 1:
                       cube.energy = e
                       cube.occupation = o
                       required_cubes.append(cube)
                       found = True

                       print("Found cube file for spin {s}, energy {e}, occupation {o}"\
                             .format(s=spin+1,e=e, o=o))
                       break
               if not found:
                   print("Missing cube file for spin {s}, energy {e}, occupation {o}"\
                           .format(s=spin+1,e=e,o=o))

# Prepare new cube file
print("\nInitializing STS cube")
stscube = cp2k.WfnCube.from_file(required_cubes[0].filename, read_data=True)

stscube.title = "STS data (z axis = energy)\n"
stscube.comment = "Range [{:4.2f} V, {:4.2f} V], de {:4.3f} V, FWHM {:4.3f} V {}\n" \
               .format(args.emin, args.emax, args.de, args.FWHM, args.bmethod)
# adjust z-dimension for energy
shape = np.array(stscube.data.shape)
shape[2] = int( (args.emax - args.emin) / args.de) + 1
stscube.data = np.zeros(shape, dtype=float)

# During export, these numbers will be
# "converted from Angstrom to Bohr"
a02A = constants.a0 / constants.Angstrom
stscube.origin[2] = args.emin * a02A
stscube.cell[2] = np.array([0,0,args.emax]) * a02A

zrange = np.linspace(args.emin, args.emax, shape[2])

if args.print_weights:
    wfile = open(args.print_weights, 'w')
    wfile.write('#index  spin   energy[eV]   weight[a.u.]\n')

# Perform STS calculation
print("\nReading data of {n} cube files".format(n=len(required_cubes)))
bar = progressbar.ProgressBar(niter=len(required_cubes))

for cube in required_cubes:

    # Reading cube files is the most time consuming part of the routine.
    # Since we need only one plane out of each cube file,
    # we save it to disk for reuse.
    planefile = "{f}.z{d}".format(f=cube.filename,d=args.height)

    # We don't want to keep all cubes in memory at once
    tmp = cp2k.WfnCube.from_cube(cube)

    plane = None
    if( os.path.isfile(planefile) ):
        plane = np.genfromtxt(planefile)

    else:
        tmp.read_cube_file(tmp.filename,read_data=True)
        if(not args.psisquared):
            tmp.data = np.square(tmp.data)

        plane = tmp.get_plane_above_atoms(args.height)
        # For STS at zero temperature, 
        # the occupation of the level in the calculation is irrelevant
        #plane = plane * tmp.occupation
        np.savetxt(planefile, plane)

    if args.print_weights:
        weight = np.sum( np.sum(plane) )
        wfile.write('{:6d} {:5d} {:10.3f} {:14.4e}\n'.\
                    format(cube.wfn, cube.spin, cube.energy, weight))

    emin = tmp.energy - eb
    emax = tmp.energy + eb
    imin = (np.abs(zrange-emin)).argmin()
    imax = (np.abs(zrange-emax)).argmin()

    for i in range(imin,imax+1):
        stscube.data[:,:,i] += plane * broadening(tmp.energy - zrange[i])

    bar.iterate()
print("\n")

if args.print_weights:
    wfile.close()

# Normalize, if asked to
if args.normalize is True:
   print("Normalizing STS data to 1")
   stscube.data /= np.sum(stscube.data)


print("\nWriting {}".format(args.outfile))
stscube.write_cube_file(args.outfile)

    
