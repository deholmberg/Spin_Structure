#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=FCup_Integrator
#SBATCH --cpus-per-task=1
#SBATCH -t4:00:00
#SBATCH --output=/farm_out/holmberg/FCup_Integrator-%j.out
#SBATCH --error=/farm_out/holmberg/FCup_Integrator-%j.err

${6}/Tree_Analysis.sh ${1} ${2} ${3} ${4} ${5}
