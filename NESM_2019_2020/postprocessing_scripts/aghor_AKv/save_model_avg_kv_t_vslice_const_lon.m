%---------------------------------
% save avg model kv at const lon
% Author: Pratik Aghor
%---------------------------------
clear all; clc; start_paths;

%---------------------------------------------------
% Specify time indices (e.g., full year)
% indxRange = 952:3877; % entire year
% indxRange=2176:2296; % Oct 1 - Oct 15, 2019
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
indxRange=3024:3152; % Jan 15 - Jan 30, 2020

Nt = length(indxRange);

% Grid and output parameters
% lon_idx_arr = [150]; % 63.25 W 
lon_idx_arr = [159]; % 63.15 W
% lon_idx_arr = [171];  % 63 W constant longitude index
lon_arr = zeros(length(lon_idx_arr), 1);  % for storing lon values
for j = 1:length(lon_idx_arr)
    lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));  % degrees east
end
lonval = lon_arr(1)

Ny = size(lat_rho_vec, 1);
N = NumLayers;  % vertical layers (Nz = N)

avg_kv_t_const_lon = zeros(N, Ny, length(lon_idx_arr));

kv_min = 1e-5; % background kv
kv_max = 5e-3; % max kv
%---------------------------------------------------
% Loop over time and accumulate AKt slices
for nt = 1:Nt
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);

    % Read AKt (vertical viscosity) and permute to [Nz+1 x Ny x Nx]
    AKt = pagetranspose(ncread(fname, 'AKt'));
    AKt = shiftdim(AKt, 2);  % (Nz+1 x Ny x Nx)

    % Interpolate AKt to rho levels
    AKt_rho = 0.5 * (AKt(1:N, :, :) + AKt(2:N+1, :, :));  % (N x Ny x Nx)
    % clip kv
    AKt_rho(AKt_rho < kv_min) = kv_min;
    AKt_rho(AKt_rho > kv_max) = kv_max;

    for j = 1:length(lon_idx_arr)
        lon_idx = lon_idx_arr(j);
        avg_kv_t_const_lon(:, :, j) = avg_kv_t_const_lon(:, :, j) + squeeze(AKt_rho(:, :, lon_idx));
    end
end

% Time average
avg_kv_t_const_lon = avg_kv_t_const_lon ./ Nt;
Z = squeeze(zr(:, :, lon_idx));
%---------------------------------------------------
% Write output to NetCDF
filename = sprintf('nesm_2019_2020_avg_kv_t_vslice_const_lon_%.2f_nt_%d_%d.nc', lonval, indxRange(1), indxRange(end));
ncid = netcdf.create(filename, 'CLOBBER');

x_len = netcdf.defDim(ncid, 'Nx', length(lon_idx_arr));
y_len = netcdf.defDim(ncid, 'Ny', Ny);
z_len = netcdf.defDim(ncid, 'Nz', N);

% Define variables
lon_idx_varid = netcdf.defVar(ncid, 'lon_idx_arr', 'double', [x_len]);
lon_varid     = netcdf.defVar(ncid, 'lon_arr',     'double', [x_len]);
kv_varid      = netcdf.defVar(ncid, 'avg_kv_t_const_lon', 'double', [z_len, y_len, x_len]);
Z_varid       = netcdf.defVar(ncid, 'Z_slice', 'double', [z_len, y_len]);
lat_vec_varid = netcdf.defVar(ncid, 'lat_rho_vec', 'double', [y_len]);

% Add metadata
netcdf.putAtt(ncid, lon_idx_varid, 'description', 'indices in lon_rho_vec where vertical slice is taken');
netcdf.putAtt(ncid, lon_varid, 'description', 'longitude values where vertical slice is taken');
netcdf.putAtt(ncid, lon_varid, 'units', 'degrees_west');

netcdf.putAtt(ncid, kv_varid, 'description', 'time-averaged vertical slice of Kv at constant lon');
netcdf.putAtt(ncid, kv_varid, 'units', 'm^2/s');
netcdf.putAtt(ncid, kv_varid, 'time_avg_indices', [indxRange(1) indxRange(end)]);

netcdf.putAtt(ncid, Z_varid, 'description', 'Z_vslice at constant lon');
netcdf.putAtt(ncid, Z_varid, 'units', 'm');

netcdf.putAtt(ncid, lat_vec_varid, 'description', 'lat_rho_vec');
netcdf.putAtt(ncid, lat_vec_varid, 'units', 'degrees N');

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid, lon_idx_varid, lon_idx_arr);
netcdf.putVar(ncid, lon_varid, lon_arr);
netcdf.putVar(ncid, kv_varid, avg_kv_t_const_lon);
netcdf.putVar(ncid, Z_varid, Z);
netcdf.putVar(ncid, lat_vec_varid, lat_rho_vec);

netcdf.close(ncid);
disp('Done saving avg_kv_t_const_lon!');
%---------------------------------------------------
