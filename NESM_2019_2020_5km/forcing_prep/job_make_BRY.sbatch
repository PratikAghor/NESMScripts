#!/bin/bash
#SBATCH -Jmake_BRY
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4G 
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH -t24:00:00
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/forcing_prep
WORKDIR=./
cd $WORKDIR
matlab -nodesktop -nosplash -nodisplay < make_BRY_hycom.m

