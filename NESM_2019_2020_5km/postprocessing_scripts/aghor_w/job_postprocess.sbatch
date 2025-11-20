#!/bin/bash
#SBATCH -Jssp_w
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=96G 
#SBATCH -t24:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/postprocessing_scripts/aghor_KmKe
WORKDIR=./
cd $WORKDIR

# 1. save rho_w annual mean
# matlab -nodesktop -nosplash -nodisplay < save_rho_w_annual_mean.m

# save_box_avg_w_t
matlab -nodesktop -nosplash -nodisplay < save_w_t_hslice.m

# save_box_avg_w_t
matlab -nodesktop -nosplash -nodisplay < save_box_avg_w_t.m

# 1. save w_t_vslice
# matlab -nodesktop -nosplash -nodisplay < save_w_t_vslice_const_lat.m
