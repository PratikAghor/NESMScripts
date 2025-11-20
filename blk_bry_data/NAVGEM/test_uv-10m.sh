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
# export  W=navgem2.0_0.176c-sec_2021_03hr  ####your file name
# export  W=navgem2.0_0.176c-sec_2020_03hr  ####your file name
export W=navgem1.4_0.281c-sec_2020_03hr
# export W=navgem2.0_0.176c-sea_2020_03hr
echo ${W}
#C
export FOR071=${W}_uv-10m.D
export CDF_FILE=${W}_uv-10m.nc
#C
/bin/rm -f $CDF_FILE
#C
./nrl2nc <<E-o-D
 &NRL2NCDF
  NFLD     = 2,
  CNAME    = 'wndewd',
             'wndnwd',
  FLG_SCL  = 1,              !0=input; 1=fit; 2=centered; 3=+ve; 4=-ve;
             1,              !0=input; 1=fit; 2=centered; 3=+ve; 4=-ve;
  T_BOUNDS = 0.0, 0.0,
             0.0, 0.0,  !monthly climatology
 &END
E-o-D
