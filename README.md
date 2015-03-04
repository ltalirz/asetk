asetk
=====

Toolkit for working with atomic and electronic structure data,
making use of the Atomistic Simulation Environment (ASE) for
atomic structures.

Formats supported:

 * IGOR Pro: Writing of .itx files
 * CP2K: Reading of energy levels in CP2K output and .MOLog files,
         reading and extrapolating of cube files.
 * Quantum ESPRESSO: Reading of some aspects from data.xml
 * Yambo: Reading of ndb.QP database as well as the yambo output and o.qp files

The toolkit is released under the MIT license.
Note that ASE is released under the GNU Lesser General Public License (LGPL).


Installation requirements
-------------------------

 * Python 2.7.5 or greater - www.python.org
 * NumPy 1.9 or greater - www.numpy.org
 * SciPy 0.14 or greater - www.scipy.org
 * matplotlib 1.4 or greater - [matplotlib.org](matplotlib.org)
 * ASE 3.8.1 or greater - [wiki.fysik.dtu.dk/ase](wiki.fysik.dtu.dk/ase)

Installation instructions
-------------------------

Let $asetk_root be the directory containing the scripts/ asetk/ subdirectories

 1. Let your python distribution know where to find the asetk module by adding 
    the $asetk_root directory to the PYTHONPATH environment variable.

    On Linux and MacOS add to your .bashrc

    ``` export PYTHONPATH=$PYTHONPATH:$asetk_root  ```
 2. Let your system know where to find the scripts by adding the scripts/
    directory to the PATH environment variable.

    On Linux and MacOS add to your .bashrc

    ``` export PATH=$PATH:$asetk_root/scripts  ```
    
