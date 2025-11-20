%--------------------------
% save daily avg speed vslice at const lon
%--------------------------
clear all; close all; clc;
start_paths; 
%--------------------------
% indxRange = 1600:1607; % July 21, 2019
% indxRange = 1600:1655; % July 21 - 27, 2019
indxRange = 2912:3023; % January 1-14, 2020
% indxRange = 3024:3079; % Jan 15 - 21, 2020
% indxRange = 3080:3135; % Jan 21 - 27, 2020

Nt = length(indxRange);

%-------------------------------------------------------
lon_idx_arr = ([159]); % 63.15 W
% lat_idx_arr = ([133]); % 38.7 N

% lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
% for i = 1:length(lat_idx_arr)
%        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
% end
% lonsec = lon_rho_vec;
% X = (repmat(lon_rho_vec, N, 1));

lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
        lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));
end
latsec = lat_rho_vec;

Y = repmat(lat_rho_vec', N, 1);

avg_speed_t_const_lon = zeros(N, Ny);
%-------------------------------------------------------
%-------------------------------------------------------
for nt = 1:Nt
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);
    hisfile = fname;

    t_arr(nt) = ncread(fname, 'time');
	
    utmp=pagetranspose(ncread(hisfile, 'u'));
    vtmp=pagetranspose(ncread(hisfile, 'v'));
    utmp = shiftdim(utmp, 2);
    vtmp = shiftdim(vtmp, 2);

    size(utmp)
    size(vtmp)

    u = u2rho_3d(utmp);
    v = v2rho_3d(vtmp);

    size(u)
    size(v)

    speed = sqrt(u.^2 + v.^2);
    
    for j = 1:1 % length(lon_idx_arr)
        lon_idx = lon_idx_arr(j);
	speed_vslice = squeeze(speed(:, :, lon_idx));
        avg_speed_t_const_lon = avg_speed_t_const_lon + speed_vslice;
    end

end

avg_speed_t_const_lon = avg_speed_t_const_lon./Nt;
Z = (squeeze(zr(:, :, lon_idx_arr(1))));

size(avg_speed_t_const_lon)
%-------------------------------------------------------


% Save avg_speed_t_const_lon vertical slice at single lon to NetCDF
% Nlon = 1 assumed
%-----------------------------------------------------------

out_file = sprintf('nesm_2019_2020_avg_speed_t_vslice_const_lon_%.2f_nt_%d_%d.nc', lon_rho_vec(lon_idx), indxRange(1), indxRange(end));
fprintf('Saving data to %s\n', out_file);

Nz = size(avg_speed_t_const_lon, 1);
Ny = size(avg_speed_t_const_lon, 2);

ncid = netcdf.create(out_file, 'NETCDF4');

% Define dimensions
dimid_z = netcdf.defDim(ncid, 'Nz', Nz);
dimid_y = netcdf.defDim(ncid, 'Ny', Ny);

% Define variables
varid_lon_idx = netcdf.defVar(ncid, 'lon_idx', 'double', []);
netcdf.putAtt(ncid, varid_lon_idx, 'description', 'Lon index in lon_rho_vec');

varid_lon = netcdf.defVar(ncid, 'lon', 'double', []);
netcdf.putAtt(ncid, varid_lon, 'description', 'Lon value (degrees W)');

varid_avg_kv = netcdf.defVar(ncid, 'avg_speed_t_const_lon', 'double', [dimid_z, dimid_y]);
netcdf.putAtt(ncid, varid_avg_kv, 'description', 'Time-averaged speed at constant lon');
netcdf.putAtt(ncid, varid_avg_kv, 'units', 'm/s');

varid_Z = netcdf.defVar(ncid, 'Z_slice', 'double', [dimid_z, dimid_y]);
netcdf.putAtt(ncid, varid_Z, 'description', 'Depth values (m) for vertical slice at constant lon');

varid_lat = netcdf.defVar(ncid, 'lat_rho_vec', 'double', dimid_y);
netcdf.putAtt(ncid, varid_lon, 'description', 'Lat vector corresponding to Ny');

netcdf.endDef(ncid);

% Write variables
netcdf.putVar(ncid, varid_lon_idx, lon_idx_arr(1));
netcdf.putVar(ncid, varid_lon, lon_arr(1));

netcdf.putVar(ncid, varid_avg_kv, avg_speed_t_const_lon);
netcdf.putVar(ncid, varid_Z, Z);

netcdf.putVar(ncid, varid_lat, lat_rho_vec);

netcdf.close(ncid);

fprintf('Done saving avg_speed_vslice_const_lon.\n');
%-----------------------------------------------------------

