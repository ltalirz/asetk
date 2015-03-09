#!/bin/bash

echo "### Reading environment variables ###"
source ../environment_variables

# where stm stm output will be stored
stm_dir='./stm'
stm_ex_dir='./stm-extrapolated'


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
        NHOMO 6
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
    ### Note ###
    # Using more accurate basis sets will extend the range of validity
    # for the tip-sample distance. Still, even using the TZV2P basis
    # set, at 6 Angstroms tip-sample distance you will observe a difference,
    # when using extrapolation.
    &KIND C
      #BASIS_SET TZV2P-MOLOPT-GTH  
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND H
      #BASIS_SET TZV2P-MOLOPT-GTH
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
EOF

echo -e "\n\n"
echo "### Running CP2K calculation ###"
echo "${cp2k_binary} < scf.inp | tee scf.out"
${para_prefix} ${cp2k_binary} -i scf.inp | tee scf.out


mkdir -p $stm_dir
cd $stm_dir

echo "### Summing cube files ###"
cp2k-sumbias.py \
  --cubes ../ANTHRACENE-WFN*cube \
  --levelsfile ../scf.out \
  --vmin -3.0 --vmax 0.0 --vstep 1.0 \
  | tee sumbias.out

echo "### Performing STM simulation ###"
stm.py \
  --stmcubes stm_*.cube \
  --heights 6.0 \
  --format 'plain' \
  --plot \
  | tee stm.out

cd ../


echo "### Extrapolating wave functions ###"
cp2k-extrapolate.py \
  --hartree ANTHRACENE-v_hartree-1_0.cube \
  --levelsfile scf.out \
  --height 2.5 \
  --extent 10 \
  ANTHRACENE*WFN*cube \
  | tee extrapolate.out

mkdir -p $stm_ex_dir
mv x.*cube $stm_ex_dir
cd $stm_ex_dir

echo "### Summing cube files ###"
cp2k-sumbias.py \
  --cubes x.ANTHRACENE-WFN*cube \
  --levelsfile ../scf.out \
  --vmin -3.0 --vmax 0.0 --vstep 1.0 \
  | tee sumbias.out

echo "### Performing STM simulation ###"
stm.py \
  --stmcubes stm_*.cube \
  --heights 6.0 \
  --format 'plain' \
  --plot \
  | tee stm.out

cd ..

echo ""
echo "### FINISHED ###"
echo "See"

echo -e "$stm_dir\t\t\t for regular STM simulation"
echo -e "$stm_ex_dir\t STM for simulation based on extrapolated wave functions"

