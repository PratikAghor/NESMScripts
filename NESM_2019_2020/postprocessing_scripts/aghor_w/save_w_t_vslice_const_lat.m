%--------------------------------------------------
% save vslice at const lat of w_t
% Author: Pratik Aghor
%--------------------------------------------------
clc;
clear all;
close all;
start_paths;
%---------------------------------------------------
%---------------------------------------------------
indxRange = 0:3877; % entire year
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
lat_idx_arr = ([112, 134]);

lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));
%---------------------------------------------------
% annual_mean_indxRange = 0:3382;
% [~, mean_Nt] = size(annual_mean_indxRange);
% annual_mean_file = strcat('nesm_2019_2020_rho_w_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

% wb = ncread(annual_mean_file, 'wb');
% rhob = ncread(annual_mean_file, 'rhob');

% size(wb)
% size(rhob)
%---------------------------------------------------
w_t_vslice_const_lat = zeros(Nt, Nz, length(lat_idx_arr), Nx);
t_arr = zeros(Nt, 1);

for nt=1:Nt
	sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
	fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
	hisfile = fname;
	
	t_arr(nt) = ncread(hisfile, 'time');

    	wtmp = zeros(N, M, L);
	wtmp = pagetranspose(ncread(hisfile, 'w'));
    	wtmp = shiftdim(wtmp, 2);
	size(wtmp)

	for i = 1:length(lat_idx_arr)
		lat_idx = lat_idx_arr(i);	
		w_t_vslice_const_lat(nt, :, i, :) = squeeze(wtmp(:, lat_idx, :));
	end
end
%--------------------------------------
% save data into a netcdf file

filename = strcat('nesm_2019_2020_w_t_vslice_const_lat_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
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

w_varid = netcdf.defVar(ncid, 'w_t_vslice_const_lat', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'timeseries of w_t at const_lat');
netcdf.putAtt(ncid, w_varid, 'units', 'm/s');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(w_t_vslice_const_lat));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, w_varid, w_t_vslice_const_lat);
% close netcdf file
netcdf.close(ncid);
%----------------------------------------
% save data into a netcdf file
disp('Done saving w_t_vslice_const_lat!')
%--------------------------------------------------------
