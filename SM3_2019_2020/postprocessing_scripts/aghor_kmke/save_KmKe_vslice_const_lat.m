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
% assumes save_KmKe_3d has saved 3d HRS, VRS, KmKe matrices
% take a vslice at const. lat values. 
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % What time indices do you need?
% indxRange = 2288:2416; % Oct. 15-30, 2019
% indxRange = 3024:3152; % Jan. 15-30, 2020

nt0=indxRange(1);
[~, Nt] = size(indxRange);

filename = strcat('sm3_2019_2020_KmKe_3d_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

KmKe3D = ncread(filename, 'KmKe');
HRS3D = ncread(filename, 'HRS');
VRS3D = ncread(filename, 'VRS');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat_idx_arr = ([111, 133]);
lat_idx_arr = ([133]);
lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
hrs_vslice = zeros(N, length(lat_arr), L);
vrs_vslice = zeros(N, length(lat_arr), L);
kmke_vslice = zeros(N, length(lat_arr), L);
%%%%%%%%%%%%%%%%
% length(lat_arr)
for i = 1:length(lat_arr) % index for layer numbers
        lat_val = lat_arr(i); % actual const latitude
        % lat_idx = find(abs(lat_rho_vec-lat_arr(i))<1e-3); % find idx of lat in lat_rho_vec
        lat_idx = lat_idx_arr(i)
        Z = (squeeze(zr(:,lat_idx,:)));
        hrs_vslice(:, i, :) = (HRS3D(:, lat_idx, :));
        vrs_vslice(:, i, :) = (VRS3D(:, lat_idx, :));
        kmke_vslice(:, i, :) = (KmKe3D(:, lat_idx, :));
end
size(hrs_vslice)
%%%%%%%%%%%%%%%%
filename = strcat('sm3_2019_2020_KmKe_vslice_const_lat_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');

% Define dimensions
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Nlat', length(lat_arr));
z_len = netcdf.defDim(ncid, 'Nz', N);
one = netcdf.defDim(ncid, 'one', 1);

% Enter define mode and define ALL variables here before endDef
hrs_varid = netcdf.defVar(ncid, 'hrs_vslice', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, hrs_varid, 'description', 'vslice at const lat horizontal Reynolds stress (HRS)');
netcdf.putAtt(ncid, hrs_varid, 'units', 'm^2s^-3');

vrs_varid = netcdf.defVar(ncid, 'vrs_vslice', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, vrs_varid, 'description', 'vslice at const lat vertical Reynolds stress (VRS)');
netcdf.putAtt(ncid, vrs_varid, 'units', 'm^2s^-3');

rs_varid = netcdf.defVar(ncid, 'kmke_vslice', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, rs_varid, 'description', 'vslice at const lat kinetic energy transfer (KmKe)');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');

lat_idx_varid = netcdf.defVar(ncid, 'lat_idx', 'double', [y_len]);
netcdf.putAtt(ncid, lat_idx_varid, 'description', 'indices in lat_rho_vec for const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_idx_varid, 'units', '--');

lat_varid = netcdf.defVar(ncid, 'lat', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');

lon_varid = netcdf.defVar(ncid, 'lon_rho_vec', 'double', [x_len]);
netcdf.putAtt(ncid, lon_varid, 'description', 'longitude vector at rho points');
netcdf.putAtt(ncid, lon_varid, 'units', 'degrees_west');

Z_varid = netcdf.defVar(ncid, 'Z_slice', 'double', [z_len x_len]); % dims order: depth x lon
netcdf.putAtt(ncid, Z_varid, 'description', 'vertical coordinate slice at constant latitude');
netcdf.putAtt(ncid, Z_varid, 'units', 'meters');

netcdf.endDef(ncid);  % ONLY ONE endDef here

% Now write variables (data arrays)
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr'); 
netcdf.putVar(ncid, lat_varid, lat_arr); 
netcdf.putVar(ncid, hrs_varid, hrs_vslice);
netcdf.putVar(ncid, vrs_varid, vrs_vslice);
netcdf.putVar(ncid, rs_varid, kmke_vslice);
netcdf.putVar(ncid, lon_varid, lon_rho_vec);

% Permute Z to match dimension order [Nz x Nx] (Z is [N x Nx])
netcdf.putVar(ncid, Z_varid, Z);

netcdf.close(ncid);

disp('Done creating HRS, VRS, KmKe vslices at const lat!')
%------------------------------------------------