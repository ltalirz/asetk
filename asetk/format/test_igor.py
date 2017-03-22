""" Tests for the IGOR format

"""
from . import igor
import unittest
import numpy as np
import numpy.testing as nt

class Wave2d(unittest.TestCase):
    """ Tests for 2d waves

    """

    def test_constructor(self):
	""" Check constructor with different parameters

	"""
        n = 100
	m = 100
        data = np.random.rand(n,m)
	xmin = 0
	xdelta = 0.1
	ymin = 0
	ydelta = 0.1

        wave1 = igor.Wave2d(
           data=data,
           xmin=xmin,
           xdelta=xdelta, 
           ymin=ymin,
           ydelta=ydelta, 
           xlabel='x [Angstroms]',
           ylabel='y [Angstroms]',
	)

        wave2 = igor.Wave2d(
           data=data,
           xmin=xmin,
           xmax=n*xdelta, 
           ymin=ymin,
           ymax=n*ydelta, 
           xlabel='x [Angstroms]',
           ylabel='y [Angstroms]',
	)

	# check that xdelta, ydelta are properly set in both cases
	for i in range(2):
	    self.assertAlmostEqual(wave1.axes[i].delta, wave2.axes[i].delta)
