#! /bin/bash
#
#set -x
set time = 1
#C
#C --- Create netCDF from a NRL precip file
#C
#C
#C --- Input.
#C
# export W=navgem2.0_0.176c-sec_2021_03hr  ####your file name
# export W=navgem2.0_0.176c-sec_2020_03hr  ####your file name
export W=navgem1.4_0.281c-sec_2020_03hr
# export W=navgem1.4_0.281c-sec_2019_03hr
#C
export FOR071=${W}_solrad.D
export CDF_FILE=${W}_solrad.nc
#C
/bin/rm -f $CDF_FILE
#C
./nrl2nc <<E-o-D
 &NRL2NCDF
  NFLD     = 1,
  CNAME    = 'shwflx',
  FLG_SCL  = 3,       
  T_BOUNDS = 0.0, 0.0,
 &END
E-o-D
