#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=800
#SBATCH --job-name=ElasticSkim
#SBATCH --cpus-per-task=1
#SBATCH -t4:30:00
#SBATCH --output=/farm_out/holmberg/ElasticSkim-%j.out
#SBATCH --error=/farm_out/holmberg/ElasticSkim-%j.err

export TMPDIR=/volatile/clas12/holmberg/junk
export TEMP=/volatile/clas12/holmberg/junk
export TMP=/volatile/clas12/holmberg/junk

echo ${1} | clas12root -l -q ../analysis_elastic.cc

#${6}Tree_Analysis.sh ${1} ${2} ${3} ${4} ${5}
