#!/bin/bash
#SBATCH -Jreynoldsstresses
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


# save speed hslice, box avg it later to find spin-up time
matlab -nodesktop -nosplash -nodisplay < save_speed_hslice.m

# save speed hslice, box avg it later to find spin-up time
matlab -nodesktop -nosplash -nodisplay < save_box_avg_speed_t.m

# from above two, determined spin-up to be ~ 120 days, so avg over 952:3877 
# save annual mean
# matlab -nodesktop -nosplash -nodisplay < save_annual_mean.m

# save eke_t_hslice at const vlevel
matlab -nodesktop -nosplash -nodisplay < save_eke_t_hslice.m
# save box_avg_eke_t at const vlevel
matlab -nodesktop -nosplash -nodisplay < save_box_avg_eke_t.m

# save eke 2d col avg
# matlab -nodesktop -nosplash -nodisplay < save_eke_t_col_avg.m
