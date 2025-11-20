% save_xyz_avg_kv_t.m
%--------------------------------------------------------------------------
% Save vertically (z=-1800 to zbottom) and box averaged Kv = <Kv> versus t
% from -1800 m to bottom, plus horizontal box averaging
% Author: Pratik Aghor
%--------------------------------------------------------------------------
clear all; close all; clc;
start_paths;  % Your setup script that defines paths, zw, etc.
%-------------------------------------------------------
indxRange = 952:3877; % entire year
% indxRange = 952:960;
Nt = length(indxRange);
%-------------------------------------------------------
lat_idx_arr = ([112]);

lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));

kv_vslice_const_lat = zeros(Nt, N, length(lat_idx_arr), Nx);
t_arr = zeros(Nt, 1);
%-------------------------------------------------------
for nt = 1:Nt
    fname = sprintf([dirHR, 'SM3_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);

    t_arr(nt) = ncread(fname, 'time');

    AKt = pagetranspose(ncread(fname, 'AKt')); % size [Nx, Ny, Nz+1]
    AKt = shiftdim(AKt, 2); % now it is (Nz+1) x Ny x Nx, like I prefer, zw also has (Nz+1) x Ny x Nx size
    
    % Interpolate AKt to rho points
    AKt_rho = 0.5 * (AKt(1:N,:,:) + AKt(2:N+1,:,:));

    if(nt == 1)
		size(AKt_rho)
    end

    for i = 1:length(lat_idx_arr)
    	lat_idx = lat_idx_arr(i);
        kv_t_vslice_const_lat(nt, :, i, :) = squeeze(AKt_rho(:, lat_idx, :));
    end

    squeeze(kv_t_vslice_const_lat(nt, :, 1, 92))
    
end
size(kv_t_vslice_const_lat)
%-------------------------------------------------------
%-------------------------------------------------------
%% Save to NetCDF

% save data into a netcdf file

filename = strcat('sm3_2019_2020_model_kv_t_vslice_const_lat_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', length(lat_idx_arr));
z_len = netcdf.defDim(ncid, 'Nz', NumLayers);
t_len = netcdf.defDim(ncid, 'Nt', Nt);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

t_varid = netcdf.defVar(ncid, 'time', 'double', [t_len]);
netcdf.putAtt(ncid, t_varid, 'description', 'time');
netcdf.putAtt(ncid, t_varid, 'units', 's');
netcdf.putAtt(ncid, t_varid, 'array dimensions', size(t_arr));

lat_idx_varid = netcdf.defVar(ncid, 'lat_idx_arr', 'double', [y_len]);
netcdf.putAtt(ncid, lat_idx_varid, 'description', 'indices in lat_rho_vec for const lat values where vertical slice is taken');
netcdf.putAtt(ncid, lat_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(lat_arr));

lat_varid = netcdf.defVar(ncid, 'lat_arr', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(lat_arr));

w_varid = netcdf.defVar(ncid, 'kv_t_vslice_const_lat', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'timeseries of kv_t at const_lat');
netcdf.putAtt(ncid, w_varid, 'units', 'm^2/s');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(kv_t_vslice_const_lat));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, w_varid, kv_t_vslice_const_lat);
% close netcdf file
netcdf.close(ncid);
%-------------------------------------------------------
