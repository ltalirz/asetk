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


Installation requirements
-------------------------

 * Python 2.7.5 or greater - www.python.org
 * NumPy 1.9 or greater - www.numpy.org
 * SciPy 0.14 or greater - www.scipy.org
 * matplotlib 1.4 or grater - matplotlib.org
 * ASE 3.8.1 or greater - wiki.fysik.dtu.dk/ase

Installation instructions
-------------------------
 
 Nothing to be done except letting your python distribution know
 where to find the asetk module.

 This can be achieved e.g. by this line in your .bashrc

 export PYTHONPATH=$PYTHONPATH:/path/to/working/copy 

