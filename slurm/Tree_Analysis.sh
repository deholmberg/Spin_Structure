#!/bin/bash

#echo ${1} ${2} ${3} ${4} | clas12root -l -q ../Get_DF_By_Sector.C

export TMPDIR=/volatile/clas12/holmberg/junk
export TEMP=/volatile/clas12/holmberg/junk
export TMP=/volatile/clas12/holmberg/junk

echo ${1} ${2} ${3} ${4} ${5} | clas12root -l -q ../GetPairTree.C
#echo ${1} ${2} ${3} | clas12root -l -q ../GetChargeFC.C
#echo ${1} ${2} ${3} ${4} ${5} | clas12root -l -q ../GetHELscalers.C

