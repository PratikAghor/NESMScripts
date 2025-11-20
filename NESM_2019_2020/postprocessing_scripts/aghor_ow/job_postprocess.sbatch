#!/bin/bash
#SBATCH -Jow_t_hslice
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=96G 
#SBATCH -t10:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/aghor_KmKe
WORKDIR=./
cd $WORKDIR

# save ow_t_hslice
matlab -nodesktop -nosplash -nodisplay < save_ow_t_hslice.m


#------------------------------------------------------------------------
# save seasonal means of vort hslice at a given vlevel
# matlab -nodesktop -nosplash -nodisplay < save_seasonal_avg_vort_hslice.m
