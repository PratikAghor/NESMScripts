#!/bin/bash
#SBATCH -Jmodel_kv_t
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=96G 
#SBATCH -t3:00:00
#SBATCH -oReport-%j.out                  # Combined output and error messages file
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
WORKDIR=./
cd $WORKDIR

# save kv_t box average at a const depth 
# matlab -nodesktop -nosplash -nodisplay < save_model_box_avg_kv_t.m
# matlab -nodesktop -nosplash -nodisplay < save_model_box_avg_kv_t_orig.m

# save xyz_avg_kv_t
# matlab -nodesktop -nosplash -nodisplay < save_model_xyz_avg_kv_t.m

# save kv_t vslice
# matlab -nodesktop -nosplash -nodisplay < save_model_kv_t_vslice_const_lat.m
# matlab -nodesktop -nosplash -nodisplay < save_model_kv_t_vslice_const_lon.m

# save_avg_kv_t_const_lon
# matlab -nodesktop -nosplash -nodisplay < save_model_avg_kv_t_vslice_const_lon.m

# save_avg_kv_t_const_lat
matlab -nodesktop -nosplash -nodisplay < save_model_avg_kv_t_vslice_const_lat.m
