#!/usr/bin/env python
import numpy as np
#import matplotlib.pyplot as plt
import argparse
import atk.format.cp2k as cp2k
import atk.util.progressbar as progressbar

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

gaussian = lambda x: 1/(args.sigma * np.sqrt(2*np.pi)) \
                     * np.exp( - x**2 / (2 * args.sigma**2) )

spectrum = cp2k.Spectrum.from_mo(args.levelsfile)
print("Read spectrum from {f} containing {s} states"\
        .format(f=args.levelsfile, s=len(spectrum.energies)))

# Reading headers of cube files
cubes = []
print("Reading headers of {n} cube files".format(n=len(args.cubes)))
bar = progressbar.ProgressBar(niter=len(args.cubes))

for fname in args.cubes:
    cubes.append( cp2k.WfnCube.from_file(fname) )
    bar.iterate()

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

# Prepare new cube file
stscube = cp2k.Cube.from_cube(required_cubes[0])

stscube.title = "STS data (z axis = energy)"
self.comment = "Range [{:4.2f} V, {4.2f} V], delta {:4.3f} V, sigma {:4.3f} V" \
               .format(args.emin, args.emax, args.de, args.sigma)
# adjust z-dimension for energy
shape = self.data.shape
shape[2] = int( (args.emax - args.emin) / args.de) + 1
self.data = np.zeros(shape)
origin[2] = emin

zrange = r_[origin[2], args.emax, shape[2]]


# Perform STS calculation
for cube in required_cubes:
    print("Processing {}".format(cube.filename))

    # Reading cube files is the most time consuming part of the routine.
    # Since we need only one plane out of each cube file,
    # we save it to disk for reuse.
    planefile = "{f}.dz{d}".format(f=cube.filename,d=args.dz)
    if( os.isfile(planefile) ):
        plane = np.genfromtxt(planefile)

    else:
        # We don't want to keep all cubes in memory at once
        tmp = cp2k.Cube.from_cube(cube)
        tmp.read_cube_file(tmp.filename,read_data=True)
        if(not args.psisquared):
            tmp.data = np.square(tmp.data)

        plane = tmp.get_plane_above_atoms(args.dz)
        np.savetxt(planefile, plane)

    emin = tmp.energy - args.sigma * args.nsigmacut
    emax = tmp.energy + args.sigma * args.nsigmacut
    imin = (np.abs(zrange-emin)).argmin()
    imax = (np.abs(zrange-emax)).argmin()

    for i in range(imin,imax+1):
        stscube.data[:,:,i] += plane * gaussian(tmp.energy - zrange[i])

print("Writing {}".format(args.outfile))
sts.cube.write_cube_file(args.outfile)


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
   
 
    
