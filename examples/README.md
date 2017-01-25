Running examples
================

 * Make sure that python is in your PATH,
   the asetk module is found in your PYTHONPATH
 * Edit the environment_variables file in the example directory
   to reflect your environment
 * execute the run.sh script in the example directory

List of examples
================

### ex01/  

    Constant-height STM simulation for anthracene.

    The tip is placed 6 Angstroms above the molecular plane.
    Results with and without wave function extrapolation are compared.

### ex02/  

    Constant-height STS simulation for anthracene.

    The tip is placed 3 Angstroms above the molecular plane.
    The STS data is saved in the cube file format, where the z dimension
    serves as the energy axis.
    For simplicity, no extrapolation is performed.

### ex03/  

    Constant-current STM simulation for anthracene.

    The tip scans over the molecule at constant local density of states.
    For simplicity, no extrapolation is performed.
