#!/bin/bash
#SBATCH -Jvort_t
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

# save vort_t_hslice
# matlab -nodesktop -nosplash -nodisplay < save_vort_t_hslice.m
# save box_avg_vort_t at const vlevel
# matlab -nodesktop -nosplash -nodisplay < save_box_avg_vort_t.m

# save vort_t_vslice_const_lon
# matlab -nodesktop -nosplash -nodisplay < save_vort_t_vslice_const_lon.m

# save vort_t_vslice_const_lat
# matlab -nodesktop -nosplash -nodisplay < save_vort_t_vslice_const_lat.m

# save uv hslice
# matlab -nodesktop -nosplash -nodisplay < save_uv_hslice_individually.m

# save vort_avg_kde_pdf
matlab -nodesktop -nosplash -nodisplay < save_vort_avg_kde_pdf.m
#------------------------------------------------------------------------
# save seasonal means of vort hslice at a given vlevel
# matlab -nodesktop -nosplash -nodisplay < save_seasonal_avg_vort_hslice.m
