#!/bin/bash


echo ${1} ${2} | clas12root -l '../RunCounts.C'

#echo ${1} ${2} ${3} ${4} ${5} | clas12root -l -q ../GetPairTree.C
#echo ${1} ${2} ${3} | clas12root -l -q ../GetChargeFC.C
#echo ${1} ${2} ${3} ${4} ${5} | clas12root -l -q ../GetHELscalers.C

