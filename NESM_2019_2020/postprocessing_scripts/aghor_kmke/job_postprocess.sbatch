#!/bin/bash
#SBATCH -Jnesm_kmke
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=192G 
#SBATCH -t96:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
WORKDIR=./
cd $WORKDIR

# save annual mean
matlab -nodesktop -nosplash -nodisplay < save_annual_mean.m

# save KmKe, HRS, VRS t at const vlevel 
matlab -nodesktop -nosplash -nodisplay < save_KmKe_t.m
