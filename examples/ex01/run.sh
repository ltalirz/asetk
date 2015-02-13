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
    BASIS_SET_FILE_NAME ${cp2k_pseudo_dir}/BASIS_MOLOPT
    POTENTIAL_FILE_NAME ${cp2k_pseudo_dir}/GTH_POTENTIALS
    &PRINT
      &E_DENSITY_CUBE
      &END E_DENSITY_CUBE
      &V_HARTREE_CUBE
      &END V_HARTREE_CUBE
      &MO_CUBES
        NHOMO 4
        NLUMO 4
        WRITE_CUBE T
      &END
    &END
    &MGRID
      CUTOFF 350
    &END
    &SCF
      EPS_SCF 1.0E-7
      &PRINT
        &RESTART
          &EACH
            QS_SCF 0
            GEO_OPT 1
          &END
          ADD_LAST NUMERIC
          FILENAME RESTART
        &END
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
      ABC [angstrom] 15 20 10
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
#&EXT_RESTART
# RESTART_FILE_NAME PROJ-1.restart
#&END
EOF

### Running CP2K calculation ###
echo "### Running CP2K calculation ###"
echo "${cp2k_binary} < scf.inp | tee scf.out"
${para_prefix} ${cp2k_binary} -i scf.inp | tee scf.out

### Extrapolating cube files ###


