#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2400
#SBATCH --job-name=DataLoop
#SBATCH --cpus-per-task=1
#SBATCH -t1:00:00
#SBATCH --output=/farm_out/holmberg/DataLoop-%j.out
#SBATCH --error=/farm_out/holmberg/DataLoop-%j.err

export TMPDIR=/volatile/clas12/holmberg/junk
export TEMP=/volatile/clas12/holmberg/junk
export TMP=/volatile/clas12/holmberg/junk

echo ${1} | clas12root -l -q ../Create_Skims/Elastic_Data_Loop.C

#${6}Tree_Analysis.sh ${1} ${2} ${3} ${4} ${5}
