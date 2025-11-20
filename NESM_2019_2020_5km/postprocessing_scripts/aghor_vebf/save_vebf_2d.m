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

% indxRange = 3080:3144;

vlevel_arr = ([-3000]);
vlevel=vlevel_arr(1);

nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);

% KmKe_3d_indxRange = indxRange; % 469:3382;
% [~, Nt] = size(KmKe_3d_indxRange);
VEBF_3d_file = strcat('nesm_2019_2020_5km_VEBF_3d_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

VEBF3D_sigma = ncread(VEBF_3d_file, 'VEBF3D');

VEBF2D_avg = zeros(Ny, Nx);
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

VEBF2D_avg = vinterp(VEBF3D_sigma,zr,vlevel);
%----------------------------------------------------------------------------------

%
% Save time average
%
% save data into a netcdf file
filename = strcat('nesm_2019_2020_5km_VEBF_2d_hslice_nt_',string(indxRange(1)), '_', ...
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

vebf_varid = netcdf.defVar(ncid, 'VEBF2D', 'double', [y_len x_len]);
netcdf.putAtt(ncid, vebf_varid, 'description', 'time avg of hslice of VEBF at a given vlevel');
netcdf.putAtt(ncid, vebf_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vebf_varid, 'array dimensions', size(VEBF2D_avg));

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
netcdf.putVar(ncid, vebf_varid, VEBF2D_avg); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating KmKe2D_avg at a given vlevel!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


