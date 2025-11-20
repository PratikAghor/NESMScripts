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
% save_KmK3_3d.m saves 3d KmKe in sigma coordinates
% interpolate to get KmKe_3d_zlevs with horizontal zlevels

clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
N = 3*N; % interpolate at 3N zlevels 

% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % What time indices do you need?
% indxRange = 2288:2416; % Oct. 15-30, 2019
% indxRange = 3024:3152; % Jan. 15-30, 2020

zlevs = linspace(-50, -5000, N);
vlevel_arr = zlevs;
% vlevel=vlevel_arr(1);

nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);

% KmKe_3d_indxRange = indxRange; % 469:3382;
% [~, Nt] = size(KmKe_3d_indxRange);
KmKe_3d_file = strcat('sm3_2019_2020_KmKe_3d_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

KmKe3D_sigma = ncread(KmKe_3d_file, 'KmKe');
HRS3D_sigma = ncread(KmKe_3d_file, 'HRS');
VRS3D_sigma = ncread(KmKe_3d_file, 'VRS');

KmKe3D_zlevs = zeros(N, Ny, Nx);
HRS3D_zlevs = zeros(N, Ny, Nx);
VRS3D_zlevs = zeros(N, Ny, Nx);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
for k = 1:N-1 
    %
    % Read the 3d matrix and do the interpolation
    %
    vlevel = vlevel_arr(k);
    KmKe3D_zlevs(k, :, :) = vinterp(KmKe3D_sigma,zr,vlevel);
    HRS3D_zlevs(k, :, :) = vinterp(HRS3D_sigma, zr, vlevel);
    VRS3D_zlevs(k, :, :) = vinterp(VRS3D_sigma, zr, vlevel);
end
%--------------------------------------------------------------------------
%
% Save KmKe average at zlevs
%
% save data into a netcdf file
filename = strcat('sm3_2019_2020_KmKe_3d_zlevs_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
z_len = netcdf.defDim(ncid, 'Nz', N);
% t_len = netcdf.defDim(ncid, 'Nt', Nt);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);


hrs_varid = netcdf.defVar(ncid, 'HRS3D_zlevs', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, hrs_varid, 'description', 'time avg of hslice of HRS at zlevels');
netcdf.putAtt(ncid, hrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, hrs_varid, 'array dimensions', size(HRS3D_zlevs));

vrs_varid = netcdf.defVar(ncid, 'VRS3D_zlevs', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, vrs_varid, 'description', 'time avg of hslice of VRS at a given zlevels');
netcdf.putAtt(ncid, vrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vrs_varid, 'array dimensions', size(VRS3D_zlevs));

rs_varid = netcdf.defVar(ncid, 'KmKe3D_zlevs', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, rs_varid, 'description', 'time avg of hslice of KmKe at a given zlevels');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, rs_varid, 'array dimensions', size(KmKe3D_zlevs));

vlevel_varid = netcdf.defVar(ncid, 'zlevs', 'double', [z_len]);
netcdf.putAtt(ncid, vlevel_varid, 'description', 'zlevels for horizontal slice (<0 means depth)');
netcdf.putAtt(ncid, vlevel_varid, 'units', 'm');
netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(zlevs'));

% z02_varid = netcdf.defVar(ncid, 'z02', 'double', one);
% netcdf.putAtt(ncid, z02_varid, 'description', 'z02 used for vertical integration');
% netcdf.putAtt(ncid, z02_varid, 'units', 'm');
% netcdf.putAtt(ncid, z02_varid, 'array dimensions', size(z02));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vlevel_varid, zlevs'); % start from 0
% netcdf.putVar(ncid, z02_varid, z02); % start from 0
netcdf.putVar(ncid, hrs_varid, HRS3D_zlevs); % start from 0
netcdf.putVar(ncid, vrs_varid, VRS3D_zlevs); % start from 0
netcdf.putVar(ncid, rs_varid, KmKe3D_zlevs); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating KmKe3D_zlevs!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


