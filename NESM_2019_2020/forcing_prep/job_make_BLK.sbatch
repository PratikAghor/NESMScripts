#!/bin/bash
#SBATCH -Jmake_BLK
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4G 
#SBATCH -t40
#SBATCH -oReport-%j.out
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
WORKDIR=./
cd $WORKDIR
matlab -nodesktop -nosplash -nodisplay < make_BLK_navgem.m

