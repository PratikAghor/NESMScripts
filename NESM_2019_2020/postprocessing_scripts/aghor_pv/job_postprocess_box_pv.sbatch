#!/bin/bash
#SBATCH -Jbox_avg_pv_t
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


# save_pv_t_hslice
# matlab -nodesktop -nosplash -nodisplay < save_pv_t_hslice.m
# save_box_avg_w_t
# matlab -nodesktop -nosplash -nodisplay < save_box_avg_pv_t.m

# 1. save pv_3d_t
matlab -nodesktop -nosplash -nodisplay < save_pv_3d_t.m
