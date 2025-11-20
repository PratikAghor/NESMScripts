#!/bin/bash
#SBATCH -Jgethycom_water_v
#SBATCH -A gts-abracco3-fy20phase1
#SBATCH -N 1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4G
#SBATCH -t240
#SBATCH -oReport-%j.out
#SBATCH --mail-user=paghor3@gatech.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# var=water_temp
# var=salinity
# var=water_u
var=water_v
# var=surf_el

# when /sur needs to be added to url
# var=u_barotropic_velocity
# var=v_barotropic_velocity
# var=ssh
# Change into working directory
WORKDIR=./${var}

cd $WORKDIR

begin_date=2020-04-26
end_date=2020-05-01

while [ "$begin_date" != "$end_date" ]
do
    echo Downloading $var on $begin_date ...
    file_name="${var}_${begin_date}.nc"
    file_url="https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0?var=$var&north=43&west=-65&east=-61&south=37&disableProjSubset=on&horizStride=1&time_start=${begin_date}T00%3A00%3A00Z&time_end=${begin_date}T21%3A00%3A00Z&timeStride=1&vertStride=1&addLatLon=true&accept=netcdf4"
   
    # for barotropic velo fields, modify url: /expt_93.0/sur
    # file_url="https://ncss.hycom.org/thredds/ncss/GLBy0.08/expt_93.0/sur?var=$var&north=43&west=-65&east=-61&south=37&disableProjSubset=on&horizStride=1&time_start=${begin_date}T00%3A00%3A00Z&time_end=${begin_date}T21%3A00%3A00Z&timeStride=1&vertStride=1&addLatLon=true&accept=netcdf4"
    wget --user-agent=" Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36" --no-check-certificate -O ${file_name} -c ${file_url}
    begin_date=`date -d "-1 days ago ${begin_date}" +%F`
done

