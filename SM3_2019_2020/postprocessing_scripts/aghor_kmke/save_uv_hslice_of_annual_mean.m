% save u v hslice at a const vlevel 
% from the 3d annual mean of uvw

%------------------------------------
start_paths;
%------------------------------------
uvw_indxRange = 952:3877;
[~, uvw_Nt] = size(uvw_indxRange);
vlevel = -4000;
uvw_file = strcat('sm3_2019_2020_uvw_annual_mean_3d_nt_', string(uvw_indxRange(1)), '_', string(uvw_indxRange(uvw_Nt)), '.nc');

ub = u2rho_3d(ncread(uvw_file, 'ub'));
vb = v2rho_3d(ncread(uvw_file, 'vb'));
u_hslice = vinterp(ub, zr, vlevel);
v_hslice = vinterp(vb, zr, vlevel);
%-------------------------------------------------------
%----------------------------------------------------%
% Save u, v hslice and grid info (zr, zw, lon, lat)
% in native CROCO layout without permuting
%----------------------------------------------------%

% Assumes variables in workspace:
% u_hslice (Ny, Nx), v_hslice (Ny, Nx)
% lat_rho, lon_rho (Ny, Nx)
% zr (N, Ny, Nx), zw (N+1, Ny, Nx)
% vlevel (scalar)

outfile = strcat('sm3_2019_2020_annual_mean_uv_hslice_vlevel_', string(vlevel), '.nc');

[N, Ny, Nx] = size(zr);
Nz_w = N + 1;

%-----------------------------------%
% Create NetCDF file
%-----------------------------------%
ncid = netcdf.create(outfile, 'NETCDF4');

% Define dimensions
dimid_x  = netcdf.defDim(ncid, 'Nx', Nx);
dimid_y  = netcdf.defDim(ncid, 'Ny', Ny);
dimid_z  = netcdf.defDim(ncid, 'N', N);
dimid_zw = netcdf.defDim(ncid, 'N_w', Nz_w);

% Define variables with correct dimension order
% (Ny, Nx) or (N, Ny, Nx)
varid_lon = netcdf.defVar(ncid, 'lon_rho', 'double', [dimid_y, dimid_x]);
varid_lat = netcdf.defVar(ncid, 'lat_rho', 'double', [dimid_y, dimid_x]);

varid_zr  = netcdf.defVar(ncid, 'zr', 'double', [dimid_z, dimid_y, dimid_x]);
varid_zw  = netcdf.defVar(ncid, 'zw', 'double', [dimid_zw, dimid_y, dimid_x]);

varid_u   = netcdf.defVar(ncid, 'u_hslice', 'double', [dimid_y, dimid_x]);
varid_v   = netcdf.defVar(ncid, 'v_hslice', 'double', [dimid_y, dimid_x]);

%-----------------------------------%
% Add variable attributes
%-----------------------------------%
netcdf.putAtt(ncid, varid_lon, 'units', 'degrees_west');
netcdf.putAtt(ncid, varid_lon, 'long_name', 'Longitude at rho-points');

netcdf.putAtt(ncid, varid_lat, 'units', 'degrees_north');
netcdf.putAtt(ncid, varid_lat, 'long_name', 'Latitude at rho-points');

% netcdf.putAtt(ncid, varid_zr, 'units', 'm');
% netcdf.putAtt(ncid, varid_zr, 'long_name', 'Depth at rho-points');
% 
% netcdf.putAtt(ncid, varid_zw, 'units', 'm');
% netcdf.putAtt(ncid, varid_zw, 'long_name', 'Depth at w-points');

netcdf.putAtt(ncid, varid_u, 'units', 'm/s');
netcdf.putAtt(ncid, varid_u, 'long_name', ['Zonal velocity at z = ', num2str(vlevel), ' m']);

netcdf.putAtt(ncid, varid_v, 'units', 'm/s');
netcdf.putAtt(ncid, varid_v, 'long_name', ['Meridional velocity at z = ', num2str(vlevel), ' m']);

%-----------------------------------%
% Exit define mode and write data
%-----------------------------------%
netcdf.endDef(ncid);

netcdf.putVar(ncid, varid_lon, lon_rho);
netcdf.putVar(ncid, varid_lat, lat_rho);
% netcdf.putVar(ncid, varid_zr,  zr);
% netcdf.putVar(ncid, varid_zw,  zw);
netcdf.putVar(ncid, varid_u,   u_hslice);
netcdf.putVar(ncid, varid_v,   v_hslice);

% Close the file
netcdf.close(ncid);

disp(['Saved NetCDF file: ', outfile]);
%-----------------------------------%