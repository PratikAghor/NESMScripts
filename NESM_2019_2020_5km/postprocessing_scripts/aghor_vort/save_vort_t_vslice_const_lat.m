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
% indxRange = 0:3382; % entire year
% indxRange=952:3877; % Jan 14 - Jan 20, 2020
indxRange=3088:3097; % Jan 23, 2020
% indxRange=3136:3145; % Jan 29, 2020
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
%------------------------------------------------------
% for vslice at const_lat
lat_idx_arr = ([23]);

lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));
%------------------------------------------------------
% for vslice at const lon
lon_idx_arr = ([35]);
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));
end
latsec = lat_rho_vec;
%---------------------------------------------------

%---------------------------------------------------
vort_t_vslice_const_lat = zeros(Nt, N, length(lat_idx_arr), Nx);
t_arr = zeros(Nt, 1);

for nt=1:Nt
	sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
	fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
	hisfile = fname;
	
	t_arr(nt) = ncread(hisfile, 'time');

    	vort_t_slices = zeros(NumLayers, Ny, Nx);

    	% Calculate vorticity for each layer
    	for vlevel = 1:NumLayers
        	[~, ~, ~, vort_t_slices(vlevel, :, :)]=get_vort(hisfile,gridfile,tindex,vlevel,coef);
    	end

	% for const lat
    	for i = 1:length(lat_arr) % index for layer numbers
        	lat_idx = lat_idx_arr(i);
        	lat_val = lat_arr(i); % actual const latitude
        	vort_t_vslice_const_lat(nt, :, i, :) = squeeze(vort_t_slices(:, lat_idx, :));
    	end

	% for const lon
    	% for j =1:length(lon_arr)
        %	lon_idx = lon_idx_arr(j);
        %	lon_val = lon_arr(j); % actual const longitude
        %	vort_t_vslice_const_lon(nt, :, :, j) = squeeze(vort_t_slices(:, :, lon_idx));
    	% end

end
%--------------------------------------
% save data into a netcdf file

filename = strcat('nesm_2019_2020_5km_vort_t_vslice_const_lat_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
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
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(lat_idx_arr));

lat_varid = netcdf.defVar(ncid, 'lat_arr', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(lat_arr));

w_varid = netcdf.defVar(ncid, 'vort_t_vslice_const_lon', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'timeseries of vort_t at const_lat');
netcdf.putAtt(ncid, w_varid, 'units', 's^{-1}');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(vort_t_vslice_const_lat));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lat_idx_varid, lat_idx_arr);
netcdf.putVar(ncid, lat_varid, lat_arr);
netcdf.putVar(ncid, w_varid, vort_t_vslice_const_lat);
% close netcdf file
netcdf.close(ncid);
%----------------------------------------
% save data into a netcdf file
disp('Done saving vort_t_vslice_const_lat!')
%--------------------------------------------------------
