#!/usr/bin/env python
import numpy as np
#import matplotlib.pyplot as plt
import argparse
import format.cp2k as cp2k

# Define command line parser
parser = argparse.ArgumentParser(
    description='Perform Scanning Tunneling Simulation.')
parser.add_argument('--version', action='version', version='%(prog)s 13.12.2013')
parser.add_argument(
    'cubes',
    nargs='+',
    metavar='cube files',
    help='The cube files to be sliced')
parser.add_argument(
    '--psisquared',
    metavar='True/False',
    default=False,
    type=bool,
    help='True, if cube files contain the density (psi squared).')
parser.add_argument(
    '--levelsfile',
    metavar='levels file',
    required=True,
    help='File containing the energy levels.')
parser.add_argument(
    '--outfile',
    metavar='output file',
    default='sts.cube',
    help='Name of cube file for STS simulation.')
parser.add_argument(
    '--emin',
    metavar='value',
    default=-3.0,
    type=float,
    help='Minimum bias for STS [V].')
parser.add_argument(
    '--emax',
    metavar='value',
    default=+3.0,
    type=float,
    help='Maximum bias for STS [V].')
parser.add_argument(
    '--de',
    metavar='value',
    default=+0.01,
    type=float,
    help='Bias step for STS [V].')
parser.add_argument(
    '--sigma',
    default=+0.2,
    type=float,
    help='sigma of Gaussian broadening [V]. FWHM = 2.355 sigma.')
parser.add_argument(
    '--nsigmacut',
    default=5,
    type=int,
    help='At a given energy E, consider levels within window [E-nsigmacut*sigma,E+nsigmacut*sigma].')
parser.add_argument(
    '--dz',
    metavar='delta z',
    default=2.5,
    type=float,
    help='The height above the topmost atom, where to extract the plane.')
parser.add_argument(
    '--vac',
    default=0.0,
    metavar='0.0 [eV]',
    help='The vacuum level. If specified, energies are given wrt the vacuum level.')
parser.add_argument(
    '--normalize',
    default=False,
    type=bool,
    help='Whether to normalize the STS intensity to 1.')
#parser.add_argument(
#    '--rep',
#    default=None,
#    nargs='+',
#    type=int, 
#    metavar='nx ny',
#    help='Number of replica along x and yi. If just one number is specified, it is taken for both x and y.')

args = parser.parse_args()


spectrum = cp2k.Spectrum.from_mo(args.levelsfile)

# Reading headers of cube files
cubes = []
for fname in args.cubes:
    cubes.append( cp2k.WfnCube.from_file(fname) )

# Connecting energies with required cube files
required_cubes = []
for spin, levels in zip(spectrum.spins, spectrum.energylevels):
    for index, e in enumerate(levels.energies):
        # If we need the cube file for this level..
        if e >= args.emin - args.sigma * args.nsigmacut and \
           e <= args.emax + args.sigma * args.nsigmacut:

               found = False
               for cube in cubes:
                   if cube.wfn == index and cube.spin == spin + 1:
                       cube.energy = e
                       required_cubes.append(cube)

                       print("Found cube file for spin {s}, energy {e}"\
                             .format(s=spin+1,e=e))
                       break
               if not found:
                   print("Missing cube file for spin {s}, energy {e}"\
                           .format(s=spin+1,e=e))

#
#
#
#if args.rep is not None:
#    if len(args.rep) == 1:
#        args.rep = [ args.rep, args.rep]
#    elif len(args.rep) !=2:
#        print('Invalid number of replicas requested')
#
#def get_energy(filename):
#    """
#    Extract energy from filename sts.-1.000.cube
#    """
#    e = re.search('(\-?\d+[\.\d]*?)\.cube', filename).group(1)
#    return float(e) - args.vac
#
#def set_lim():
#    ring_center = [33.6,19.1]
#    ring_diameter = 38.6
#    plt.xlim([ ring_center[0]-ring_diameter/2, ring_center[0]+ring_diameter/2])
#    plt.ylim([ ring_center[1]-ring_diameter/2, ring_center[1]+ring_diameter/2])
#
#files = args.cubes
#fig = plt.figure()
#
#for f in files:
#    plt.clf()
#    resampled, extent = sts.get_plane(f, deltaz=args.dz, rep=args.rep, nsamples=args.samples)
#
#    print('max {m}, min {min}'.format(m=np.max(resampled), min=np.min(resampled)))
#    cax = plt.imshow(resampled, extent=extent, vmax = args.vmax, cmap='gray')
#    plt.xlabel('x [$\AA$]')
#    plt.ylabel('y [$\AA$]')
#    plt.title('E={e:4.2f} eV'.format(e=get_energy(f)))
#
#    cbar = fig.colorbar(cax, format='%.2e')
#    cbar.set_label('$\\rho(E)$')
#    set_lim()
#   
#    #io.write('atoms.png',atoms)
#    #model = plt.imread('atoms.png')
#    #plt.imshow(model, extent=extent)
#    
#    outname='slice_{d:.2f}.png'.format(d=get_energy(f))
#    plt.savefig(outname, dpi=300)
#    print('Done with {f}'.format(f=f))
#  
##plt.show()
   
 
    
