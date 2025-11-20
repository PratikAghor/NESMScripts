% 
% Compute the kinetic energy transfer: horizontal Reynolds stress (HRS)
% - Pratik Aghor, modified from get_KmKe.m
% HRS= -[ <up up> dubar/dx + <up vp>  dubar/dy  ....
%          <vp up> dvbar/dx + <vp vp>  dvbar/dy ]
%     =  -<up(up.grad(ubar))+vp(up.grad(vbar))>
%
% Advection operators are used. Much easier in sigma coordinates.
%
% This is the method used in
% Djakouré, S., P. Penven, B. Bourlès, J. Veitch and V. Koné, 
% Coastally trapped eddies in the north of the Gulf of Guinea, 2014, J. Geophys. Res,
% DOI: 10.1002/2014JC010243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Assume mean and KmKe_3d with respect to that mean has been calculated and saved
% only take a section of the saved KmKe_3d, vertically integrate to take mean b/w two verical levels
% and save KmKe2D or take an hslice at a given vlevel (depth)
 
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % What time indices do you need?
% indxRange = 2288:2416; % Oct. 15-30, 2019
% indxRange = 3024:3152; % Jan. 15-30, 2020

vlevel_arr = ([-3000]);
vlevel=vlevel_arr(1);

nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);

% [~, Nt] = size(KmKe_3d_indxRange);
KmKe_3d_file = strcat('nesm_2019_2020_5km_KmKe_3d_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

KmKe3D_avg = ncread(KmKe_3d_file, 'KmKe');
HRS3D_avg = ncread(KmKe_3d_file, 'HRS');
VRS3D_avg = ncread(KmKe_3d_file, 'VRS');

KmKe2D_avg = zeros(Ny, Nx);
HRS2D_avg = zeros(Ny, Nx);
VRS2D_avg = zeros(Ny, Nx);
%----------------------------------------------------------------------------------

%
% Verticaly integrate
%
time=0;
fname = sprintf([dirHR, 'nesm_2019_2020_5km_avg.%05d.nc'], time);
hisfile = fname;
zeta = pagetranspose(ncread(hisfile, 'zeta'));

zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
% mask=maskr;
% mask(mask==0)=NaN;
% [HRS2D,h0]=vintegr2(HRS,zw,zr,z01,z02);
%----------------------------------------------------------------------------------
% take a horizontal slice at a given vlevel, modified from get_hslice.m
%
% Read the 3d matrix and do the interpolation
%
KmKe_sigma=KmKe3D_avg;
KmKe2D_avg = vinterp(KmKe_sigma,zr,vlevel);
HRS2D_avg = vinterp(HRS3D_avg, zr, vlevel);
VRS2D_avg = vinterp(VRS3D_avg, zr, vlevel);
%----------------------------------------------------------------------------------

%
% Save time average
%
% save data into a netcdf file
filename = strcat('nesm_2019_2020_5km_KmKe_2d_hslice_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '_vlevel_', string(vlevel), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% t_len = netcdf.defDim(ncid, 'Nt', Nt);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);


hrs_varid = netcdf.defVar(ncid, 'HRS2D', 'double', [y_len x_len]);
netcdf.putAtt(ncid, hrs_varid, 'description', 'time avg of hslice of HRS at a given vlevel');
netcdf.putAtt(ncid, hrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, hrs_varid, 'array dimensions', size(HRS2D_avg));

vrs_varid = netcdf.defVar(ncid, 'VRS2D', 'double', [y_len x_len]);
netcdf.putAtt(ncid, vrs_varid, 'description', 'time avg of hslice of VRS at a given vlevel');
netcdf.putAtt(ncid, vrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vrs_varid, 'array dimensions', size(VRS2D_avg));

rs_varid = netcdf.defVar(ncid, 'KmKe2D', 'double', [y_len x_len]);
netcdf.putAtt(ncid, rs_varid, 'description', 'time avg of hslice of KmKe at a given vlevel');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, rs_varid, 'array dimensions', size(KmKe2D_avg));

vlevel_varid = netcdf.defVar(ncid, 'vlevel', 'double', one);
netcdf.putAtt(ncid, vlevel_varid, 'description', 'vlevel for horizontal slice (<0 means depth)');
netcdf.putAtt(ncid, vlevel_varid, 'units', 'm');
netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(vlevel_arr));

% z02_varid = netcdf.defVar(ncid, 'z02', 'double', one);
% netcdf.putAtt(ncid, z02_varid, 'description', 'z02 used for vertical integration');
% netcdf.putAtt(ncid, z02_varid, 'units', 'm');
% netcdf.putAtt(ncid, z02_varid, 'array dimensions', size(z02));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vlevel_varid, vlevel); % start from 0
% netcdf.putVar(ncid, z02_varid, z02); % start from 0
netcdf.putVar(ncid, hrs_varid, HRS2D_avg); % start from 0
netcdf.putVar(ncid, vrs_varid, VRS2D_avg); % start from 0
netcdf.putVar(ncid, rs_varid, KmKe2D_avg); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating KmKe2D_avg at a given vlevel!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


