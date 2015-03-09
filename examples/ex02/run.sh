#!/bin/bash

echo "### Reading environment variables ###"
source ../environment_variables

cat > scf.inp <<EOF
&GLOBAL
  PROJECT ANTHRACENE
  PRINT_LEVEL LOW
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME ${cp2k_libdir}/BASIS_MOLOPT
    POTENTIAL_FILE_NAME ${cp2k_libdir}/GTH_POTENTIALS
    &PRINT
      &E_DENSITY_CUBE
      &END E_DENSITY_CUBE
      &V_HARTREE_CUBE
      &END V_HARTREE_CUBE
      &MO_CUBES
        NHOMO 7
        NLUMO 3
        WRITE_CUBE T
      &END
    &END
    &MGRID
      CUTOFF 300
    &END
    &SCF
      EPS_SCF 1.0E-7
      &PRINT
        &RESTART_HISTORY OFF
        &END
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 15 20 15
    &END
    &TOPOLOGY
     COORD_FILE_NAME ./anthracene.xyz
     COORDINATE xyz
    &END
    &KIND C
      BASIS_SET TZV2P-MOLOPT-GTH  
      #BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      BASIS_SET TZV2P-MOLOPT-GTH
      #BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
EOF

echo -e "\n\n"
echo "### Running CP2K calculation ###"
echo "${cp2k_binary} < scf.inp | tee scf.out"
#${para_prefix} ${cp2k_binary} -i scf.inp | tee scf.out


echo "### Performing STS simulation ###"
sts.py \
  --cubes ANTHRACENE-WFN*cube \
  --levelsfile scf.out \
  --vmin -3.0 --vmax 1.0 --vstep 0.02\
  --FWHM 0.1\
  --bmethod Gaussian\
  --height 3.0\
  | tee sts.out

echo -e "\n"
echo "### Extracting plane for plotting ###"
cube-ex-plane.py sts.cube x 16 

echo "python ./plot.py"
python ./plot.py

echo -e "\n"
echo "### FINISHED ###"
echo "See"
echo -e "sts.cube\t\t for STS simulation (z-dimension = voltage)"
echo -e "sts.cube.plane16.dat\t for plane extracted from cube file (x index 16)."
echo -e "sts.cube.plane16.png\t for plot of extracted plane"

