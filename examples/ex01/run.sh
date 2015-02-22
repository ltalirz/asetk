#!/bin/bash

echo "### Reading environment variables ###"
source ./environment_variables

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
        NHOMO 5
        NLUMO 5
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
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      BASIS_SET TZV2P-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
EOF

echo "### Running CP2K calculation ###"
echo "${cp2k_binary} < scf.inp | tee scf.out"
${para_prefix} ${cp2k_binary} -i scf.inp | tee scf.out

echo "### Summing cube files ###"
cp2k-sumbias.py \
  --cubes ANTHRACENE-WFN*cube \
  --levelsfile scf.out \
  --vmin -3.0 --vmax 1.0 --vstep 1.0 \
  | tee sumbias.out

echo "### Performing STM simulation ###"
stm.py \
  --stmcubes stm_*.cube \
  --heights 6.0 \
  --format 'plain' \
  --plot \
  | tee stm.out


echo "### Extrapolating wave functions ###"
cp2k-extrapolate.py \
  --hartree ANTHRACENE-v_hartree-1_0.cube \
  --levelsfile scf.out \
  --height 2.5 \
  --extent 10 \
  ANTHRACENE*WFN*cube \
  | tee extrapolate.out

mkdir -p extrapolated
mv x.*cube extrapolated/
cd extrapolated

echo "### Summing cube files ###"
cp2k-sumbias.py \
  --cubes x.ANTHRACENE-WFN*cube \
  --levelsfile ../scf.out \
  --vmin -3.0 --vmax 1.0 --vstep 1.0 \
  | tee sumbias.out

echo "### Performing STM simulation ###"
stm.py \
  --stmcubes stm_*.cube \
  --heights 6.0 \
  --format 'plain' \
  --plot \
  | tee stm.out

cd ..

echo "### FINISHED ###"
echo "Find regular STM simulation in current working directory"
echo "and extrapolated results in extrapolated/"

