#!/bin/bash
#SBATCH -Jmerge_hycom_data
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4G
#SBATCH -oReport_%j.out
#SBATCH -t40
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Change into working directory
# WORKDIR=/storage/home/hcoda1/9/paghor3/scratch/work/SEAMOUNT_2019_2020/data/merged
WORKDIR=./
cd $WORKDIR

# HYCOMDIR=/storage/home/hcoda1/9/paghor3/scratch/work/ATLANTIS_2019_2020/data/HYCOM
HYCOMDIR=../HYCOM

begin_date=20190101
end_date=20200501

# varArray=("salinity" "water_temp" "water_u" "water_v" "surf_el")
varArray=("water_v")
for i in ${!varArray[@]}; do
	var=${varArray[$i]}
	echo "merging ${var} data" \n 
	cdo mergetime ${HYCOMDIR}/${var}/*.nc ${var}_${begin_date}_${end_date}.nc
done
