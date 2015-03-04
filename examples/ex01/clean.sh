#!/bin/bash
echo "### Cleaning all output data###"

. environment_variables

rm *.cube *.png *.dat *.out *.inp 
rm ANTHRACENE* 
rm -r $stm_dir $stm_ex_dir
