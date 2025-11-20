% build frc file for GOM_hycom run
clear;
Yorig = 2018;
dt = 1; % To build 3hr forcing: for 3hr data, use dt = 1; for 1hr data, use dt = 3;
tnum_start = datenum('2019-01-01');
tnum_end = datenum('2020-05-01');
L=299; M=277;
% L=60; M=56;
Lp=L+1;
Mp=M+1;


dir_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/NESM_2019_2020/'
blk_bry_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/blk_bry_data/'

grdname=[dir_path, 'forcing/','NESM_grd.nc'];

wndfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sec_2019_03hr_uv-10m.nc'];
wndfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sec_2020_03hr_uv-10m.nc'];

swradfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sec_2019_03hr_solrad.nc'];
swradfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sec_2020_03hr_solrad.nc'];

dlwsfcfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2019_03hr_dlwsfc.nc'];
dlwsfcfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2020_03hr_dlwsfc.nc'];

precfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-std_2019_03hr_precip.nc'];
precfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-std_2020_03hr_precip.nc'];

tairfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2019_03hr_temp2m.nc'];
tairfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2020_03hr_temp2m.nc'];

spchumfile1 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2019_03hr_spchum.nc'];
spchumfile2 = [blk_bry_path, 'NAVGEM/navgem1.4_0.281c-sea_2020_03hr_spchum.nc'];

blkout = [dir_path, 'forcing/NESM_blk_20190101_20200501.nc'];

ncid=netcdf.create(blkout,'NETCDF4');
xi_u=netcdf.defDim(ncid,'xi_u',L);
xi_v=netcdf.defDim(ncid,'xi_v',Lp);
eta_u=netcdf.defDim(ncid,'eta_u',Mp);
eta_v=netcdf.defDim(ncid,'eta_v',M);
xi_rho=netcdf.defDim(ncid,'xi_rho',Lp);
eta_rho=netcdf.defDim(ncid,'eta_rho',Mp);
netcdf.endDef(ncid);
netcdf.close(ncid)

% read model grid
lonr=ncread(grdname,'lon_rho'); latr=ncread(grdname,'lat_rho');
lonu=ncread(grdname,'lon_u'); latu=ncread(grdname,'lat_u');
lonv=ncread(grdname,'lon_v'); latv=ncread(grdname,'lat_v');
maskr=ncread(grdname,'mask_rho'); masku=ncread(grdname,'mask_u');
maskv=ncread(grdname,'mask_v');
% build wind stress forcing
disp('Working on wind stress ...')
Code_GOM_navgem_wnd

% build solar shortwave radiation
disp('Working on shortwave radiation...')
Code_GOM_navgem_radsw

% build longwave radiation
disp('Working on longwave radiation ...')
Code_GOM_navgem_radlw_in

% build precipitation rate
disp('Working on precipitation rate...')
Code_GOM_navgem_prate

% build surface air temperature
disp('Working on surface air temperature ...')
Code_GOM_navgem_tair

% build relative humidity
disp('Working on relative humidity ...')
Code_GOM_navgem_rhum
%make_rivers
