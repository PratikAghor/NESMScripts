% read the coordinate
clear;
ncfile = '../grid_1km_z100.nc';
grd.lat_rho = ncread(ncfile,'lat_rho');
grd.lon_rho = ncread(ncfile,'lon_rho');
grd.mask_rho = ncread(ncfile,'mask_rho');
grd.lat_psi = ncread(ncfile,'lat_psi');
grd.lon_psi = ncread(ncfile,'lon_psi');
grd.mask_psi = ncread(ncfile,'mask_psi');
grd.lat_u = ncread(ncfile,'lat_u');
grd.lon_u = ncread(ncfile,'lon_u');
grd.mask_u = ncread(ncfile,'mask_u');
grd.lat_v = ncread(ncfile,'lat_v');
grd.lon_v = ncread(ncfile,'lon_v');
grd.mask_v = ncread(ncfile,'mask_v');
grd.angle = ncread(ncfile,'angle');
grd.h = ncread(ncfile,'h');
grd.pm = ncread(ncfile,'pm');
grd.pn = ncread(ncfile,'pn');
grd.f = ncread(ncfile,'f');

save  grd_NESM_1km_z100.mat grd
