%-----------------------------------
% save model avg kv const lat
% Note: can be generalized, but for now,
% works only for one const lat slice.
% - Pratik Aghor
%-----------------------------------
clear all; close all; clc;
start_paths;  % Your setup script that defines paths, zw, etc.
%-------------------------------------------------------
indxRange = 952:3877; % entire year
% indxRange = 952:960;
Nt = length(indxRange);
%-------------------------------------------------------
% lat_idx_arr = ([111]); % 38.5 N
lat_idx_arr = ([133]); % 38.7 N

lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));

avg_kv_t_const_lat = zeros(N, Nx);
t_arr = zeros(Nt, 1);

kv_min = 1e-5; % background kv
kv_max = 5e-3; % max kv
%-------------------------------------------------------
for nt = 1:Nt
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);

    t_arr(nt) = ncread(fname, 'time');

    AKt = pagetranspose(ncread(fname, 'AKt')); % size [Nx, Ny, Nz+1]
    AKt = shiftdim(AKt, 2); % now it is (Nz+1) x Ny x Nx, like I prefer, zw also has (Nz+1) x Ny x Nx size

    % Interpolate AKt to rho points
    AKt_rho = 0.5 * (AKt(1:N,:,:) + AKt(2:N+1,:,:));
    
    % clip kv
    AKt_rho(AKt_rho < kv_min) = kv_min;
    AKt_rho(AKt_rho > kv_max) = kv_max;

    if(nt == 1)
                size(AKt_rho)
    end

    for i = 1:1 % length(lat_idx_arr)
        lat_idx = lat_idx_arr(i);
        avg_kv_t_const_lat = avg_kv_t_const_lat + squeeze(AKt_rho(:, lat_idx, :));
    end

end
avg_kv_t_const_lat = avg_kv_t_const_lat./Nt;
Z = (squeeze(zr(:, lat_idx_arr(1), :)));

size(avg_kv_t_const_lat)
%-------------------------------------------------------
%-----------------------------------------------------------
% Save avg_kv_t_const_lat vertical slice at single latitude to NetCDF
% Nlat = 1 assumed
%-----------------------------------------------------------

out_file = sprintf('nesm_2019_2020_avg_kv_t_vslice_const_lat_nt_%d_%d.nc', indxRange(1), indxRange(end));
fprintf('Saving data to %s\n', out_file);

Nz = size(avg_kv_t_const_lat, 1);
Nx = size(avg_kv_t_const_lat, 2);

ncid = netcdf.create(out_file, 'NETCDF4');

% Define dimensions
dimid_z = netcdf.defDim(ncid, 'Nz', Nz);
dimid_x = netcdf.defDim(ncid, 'Nx', Nx);

% Define variables
varid_lat_idx = netcdf.defVar(ncid, 'lat_idx', 'double', []);
netcdf.putAtt(ncid, varid_lat_idx, 'description', 'Latitude index in lat_rho_vec');

varid_lat = netcdf.defVar(ncid, 'lat', 'double', []);
netcdf.putAtt(ncid, varid_lat, 'description', 'Latitude value (degrees N)');

varid_avg_kv = netcdf.defVar(ncid, 'avg_kv_t_const_lat', 'double', [dimid_z, dimid_x]);
netcdf.putAtt(ncid, varid_avg_kv, 'description', 'Time-averaged vertical diffusivity Kv at constant latitude');
netcdf.putAtt(ncid, varid_avg_kv, 'units', 'm^2/s');

varid_Z = netcdf.defVar(ncid, 'Z_slice', 'double', [dimid_z, dimid_x]);
netcdf.putAtt(ncid, varid_Z, 'description', 'Depth values (m) for vertical slice at constant latitude');

varid_lon = netcdf.defVar(ncid, 'lon_rho_vec', 'double', dimid_x);
netcdf.putAtt(ncid, varid_lon, 'description', 'Longitude vector corresponding to Nx');

netcdf.endDef(ncid);

% Write variables
netcdf.putVar(ncid, varid_lat_idx, lat_idx_arr(1));
netcdf.putVar(ncid, varid_lat, lat_arr(1));

netcdf.putVar(ncid, varid_avg_kv, avg_kv_t_const_lat);
netcdf.putVar(ncid, varid_Z, Z);

netcdf.putVar(ncid, varid_lon, lon_rho_vec);

netcdf.close(ncid);

fprintf('Done saving avg_kv_vslice_const_lat.\n');
%-----------------------------------------------------------
