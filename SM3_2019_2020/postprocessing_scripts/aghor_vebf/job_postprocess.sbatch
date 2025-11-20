#!/bin/bash
#SBATCH -Jsm3_vebf
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=192G 
#SBATCH -t96:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/aghor_KmKe
WORKDIR=./
cd $WORKDIR

# save annual mean
# matlab -nodesktop -nosplash -nodisplay < save_rho_w_annual_mean.m

# save vebf t at const vlevel
matlab -nodesktop -nosplash -nodisplay < save_vebf_t.m

