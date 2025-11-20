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
% case_type ='instant'; % 'instant' or 'avg'
case_type ='avg'; % 'instant' or 'avg'

if(strcmp(case_type, 'instant'))
    nt = 3040;
elseif(strcmp(case_type, 'avg'))
    indxRange = 952:3877; % What time indices do you need?
    % indxRange = 2288:2416; % Oct. 15-30, 2019
    % indxRange = 3024:3152; % Jan. 15-30, 2020

    nt0=indxRange(1);
    [~, Nt] = size(indxRange);
end


if(strcmp(case_type, 'instant'))
    filename = strcat('nesm_2019_2020_instant_VEBF_3d_nt_',string(nt), '.nc');

elseif(strcmp(case_type, 'avg'))
    filename = strcat('nesm_2019_2020_VEBF_3d_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
end

VEBF3D = ncread(filename, 'VEBF3D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat_idx_arr = ([111, 133]);
lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
vebf_vslice = zeros(N, length(lat_arr), L);
%%%%%%%%%%%%%%%%
% length(lat_arr)
for i = 1:length(lat_arr) % index for layer numbers
        lat_val = lat_arr(i); % actual const latitude
        % lat_idx = find(abs(lat_rho_vec-lat_arr(i))<1e-3); % find idx of lat in lat_rho_vec
        lat_idx = lat_idx_arr(i)
        % Normalize with coriolis freq.
        vebf_vslice(:, i, :) = (VEBF3D(:, lat_idx, :));
end
size(vebf_vslice)
%%%%%%%%%%%%%%%%
if(strcmp(case_type, 'instant'))
    filename = strcat('nesm_2019_2020_instant_VEBF_vslice_const_lat_nt_',string(nt), '.nc');

elseif(strcmp(case_type, 'avg'))
    filename = strcat('nesm_2019_2020_VEBF_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
        string(indxRange(Nt)), '.nc');
end
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Nlat', length(lat_arr));
z_len = netcdf.defDim(ncid, 'Nz', N);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

rs_varid = netcdf.defVar(ncid, 'vebf_vslice', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, rs_varid, 'description', 'vslice at const lat VEBF');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, rs_varid, 'array dimensions', size(vebf_vslice));

lat_idx_varid = netcdf.defVar(ncid, 'lat_idx', 'double', [y_len]);
netcdf.putAtt(ncid, lat_idx_varid, 'description', 'indices in lat_rho_vec for const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(lat_idx_arr'));

lat_varid = netcdf.defVar(ncid, 'lat', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical epv slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(lat_arr));% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr'); % start from 0
netcdf.putVar(ncid, lat_varid, lat_arr); % start from 0
netcdf.putVar(ncid, rs_varid, vebf_vslice); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating VEBF vslices at const lat!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
