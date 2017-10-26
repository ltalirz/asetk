asetk
=====

Toolkit for working with atomic and electronic structure data,
built on top of the Atomistic Simulation Environment (ASE) for
atomic structures.

Formats supported:

 * CP2K: STM simulation, STS simulation, wave function extrapolation,
   reading of energy levels from CP2K output, .MOLog and .pdos files
 * IGOR Pro: Reading and writing of .itx files
 * Quantum ESPRESSO: Reading of some aspects from data.xml,
      convert intermediate format to .cube
 * Yambo: Reading of ndb.QP, yambo output and o.qp files
 * BerkeleyGW: Reading parts of eps0mat.h5
 * Gaussian Cube Format: Reading, writing, slicing, averaging, rolling, plotting

See the `scripts/` subdirectory for all stand-alone command line scripts
(and run `<script>.py -h` to find out what it does).


Installation instructions
-------------------------

```bash
git clone https://github.com/ltalirz/asetk
cd asetk
pip install -e .
```

License information
-------------------

The toolkit is released under the MIT license.
Note that ASE is released under the GNU Lesser General Public License (LGPL).

If your scientific publication has benefited from the use of asetk,
please consider an acknowledgement through the following citation

*Talirz, L., asetk; https://github.com/ltalirz/asetk, [![DOI](https://zenodo.org/badge/15176282.svg)](https://zenodo.org/badge/latestdoi/15176282)*

Contact information
-------------------

Contact [leopold.talirz@gmail.com](mailto:leopold.talirz@gmail.com) for any
questions.
