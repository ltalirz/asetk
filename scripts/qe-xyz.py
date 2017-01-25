#!/usr/bin/env python
import asetk.format.qe as qe
import sys

if len(sys.argv) >= 2:
    prefix = sys.argv[1]
else:
    raise IOError("Specify PREFIX of .save directory")

at = qe.Atoms.from_save(prefix)
at.write(prefix + ".xyz")
