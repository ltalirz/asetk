#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

fname = 'sts.cube.plane16'
data = np.genfromtxt(fname+'.dat')

# for plotting with matplotlib
data = data.swapaxes(0,1)[::-1]

fig = plt.Figure(figsize=(4,4))
#ax = fig.add_subplot(111, aspect=0.5)
#fig.add_subplot(111,aspect='equal')
plt.xlabel('y [$\AA$]')
plt.ylabel('U [V]')

extent = [0,20,-3,1]

cax = plt.imshow(data, extent=extent, cmap='gray',
             aspect=4.0)

#cbar = fig.colorbar(cax, format='%.2e')
#cbar.set_label('LDOS')

outfile = fname+'.png'
print("Plotting into {}".format(outfile))
plt.savefig(outfile, dpi=200)
