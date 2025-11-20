#! /bin/bash
#
#set -x
# usage bash test_temp2m.sh 2021
# set time=1
#C
#C --- Create netCDF from a NRL precip file
#C
#C
#C --- Input.
#C
# set year=1 # first trailing arg to sh file is year
# echo ${year}
# export  W=navgem2.0_0.176c-sea_2021_03hr  ####your file name
# export  W=navgem2.0_0.176c-sea_2020_03hr  ####your file name
export W=navgem1.4_0.281c-sea_2020_03hr
# export W=navgem1.4_0.281c-sea_2019_03hr
#C
export FOR071=${W}_temp2m.D
export CDF_FILE=${W}_temp2m.nc

echo ${FOR071}
#C
/bin/rm -f $CDF_FILE
#C
./nrl2nc <<E-o-D
 &NRL2NCDF
  NFLD     = 1,
  CNAME    = 'airtmp',
  FLG_SCL  = 1,       
  T_BOUNDS = 0.0, 0.0,
 &END
E-o-D
