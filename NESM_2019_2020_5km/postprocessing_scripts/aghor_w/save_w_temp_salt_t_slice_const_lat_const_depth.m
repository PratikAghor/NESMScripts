%--------------------------------------------------
% save slice at const lat, const depth of w_t
% Author: Pratik Aghor
%--------------------------------------------------
clc;
clear all;
close all;
start_paths;
%---------------------------------------------------
%---------------------------------------------------
indxRange = 952:3877; % entire year
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% lat_idx = 111; % 38.5 N for 1km res
% lat_idx = 23; % 38.5 N for 5km res
% lat_idx = 133; % 38.7 N for 1km res
lat_idx = 27; % 38.7 N for 5km res
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));

vlevel = -4000;
disp(['Selected vertical level: ' num2str(vlevel) ' m']);
%---------------------------------------------------

%---------------------------------------------------
w_t_slice_const_lat_const_depth = zeros(Nt, Nx);
temp_t_slice_const_lat_const_depth = zeros(Nt, Nx);
salt_t_slice_const_lat_const_depth = zeros(Nt, Nx);
dTdz_t_slice_const_lat_const_depth = zeros(Nt, Nx);
t_arr = zeros(Nt, 1);

for nt=1:Nt
	sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
	fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
	hisfile = fname;
	
	t_arr(nt) = ncread(hisfile, 'time');

    	wtmp = zeros(N, M, L);
	wtmp = pagetranspose(ncread(hisfile, 'w'));
	temp = pagetranspose(ncread(hisfile, 'temp'));
        salt = pagetranspose(ncread(hisfile, 'salt'));

	wtmp = shiftdim(wtmp, 2);
        temp = shiftdim(temp, 2);
        salt = shiftdim(salt, 2);
    
	size(wtmp)

	[dTdz_3d_zw, ~] = vert_grad(temp, zr);
	dTdz_3d = zeros(N, M, L);
	dTdz_3d(1:N-1, :, :) = dTdz_3d_zw;
        dTdz_3d(N, :, :) = dTdz_3d(N-1, :, :); % just copy the second last row, inconsequential here	
	size(dTdz_3d)

	
	% get 2d hslice
	w_2d_hslice = vinterp(wtmp, zr, vlevel);
	temp_2d_hslice = vinterp(temp, zr, vlevel);
	salt_2d_hslice = vinterp(salt, zr, vlevel);
	dTdz_2d_hslice = vinterp(dTdz_3d, zr, vlevel);

	% get 1d slice at lat_idx
	w_t_slice_const_lat_const_depth(nt, :) = w_2d_hslice(lat_idx, :);
	temp_t_slice_const_lat_const_depth(nt, :) = temp_2d_hslice(lat_idx, :);
	salt_t_slice_const_lat_const_depth(nt, :) = salt_2d_hslice(lat_idx, :);
	dTdz_t_slice_const_lat_const_depth(nt, :) = dTdz_2d_hslice(lat_idx, :);
	% ssp_t_slice_const_lat_const_depth(nt, :) = ssp_2d_hslice(lat_idx, :);

end
%--------------------------------------
% save data into a netcdf file
lat_str = sprintf('%.2f', lat_rho_vec(lat_idx));

filename = strcat('nesm_2019_2020_5km_w_temp_salt_slice_const_lat_', lat_str, '_const_z_', num2str(vlevel), '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', 1);
z_len = netcdf.defDim(ncid, 'Nz', N);
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
netcdf.putAtt(ncid, lat_idx_varid, 'array dimensions', size(lat_idx));

lat_varid = netcdf.defVar(ncid, 'lat_arr', 'double', [y_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'const lat values where vertical slice is taken');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(lat_idx));

lon_varid = netcdf.defVar(ncid, 'lon_rho_vec', 'double', [x_len]);
netcdf.putAtt(ncid, lat_varid, 'description', 'lon_rho_vec');
netcdf.putAtt(ncid, lat_varid, 'units', '--');
netcdf.putAtt(ncid, lat_varid, 'array dimensions', size(lon_rho_vec));

w_varid = netcdf.defVar(ncid, 'w_t_vslice_const_lat_const_z', 'double', [t_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'w_t at const_lat, const depth');
netcdf.putAtt(ncid, w_varid, 'units', 'm/s');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(w_t_slice_const_lat_const_depth));

temp_varid = netcdf.defVar(ncid, 'temp_t_vslice_const_lat_const_z', 'double', [t_len x_len]);
netcdf.putAtt(ncid, temp_varid, 'description', 'temp_t at const_lat, const depth');
netcdf.putAtt(ncid, temp_varid, 'units', 'deg C');
netcdf.putAtt(ncid, temp_varid, 'array dimensions', size(temp_t_slice_const_lat_const_depth));

salt_varid = netcdf.defVar(ncid, 'salt_t_vslice_const_lat_const_z', 'double', [t_len x_len]);
netcdf.putAtt(ncid, salt_varid, 'description', 'salt_t at const_lat, const depth');
netcdf.putAtt(ncid, salt_varid, 'units', 'PSU');
netcdf.putAtt(ncid, salt_varid, 'array dimensions', size(salt_t_slice_const_lat_const_depth));

dTdz_varid = netcdf.defVar(ncid, 'dTdz_t_vslice_const_lat_const_z', 'double', [t_len x_len]);
netcdf.putAtt(ncid, dTdz_varid, 'description', 'dTdz_t at const_lat, const depth');
netcdf.putAtt(ncid, dTdz_varid, 'units', 'deg C/m');
netcdf.putAtt(ncid, dTdz_varid, 'array dimensions', size(dTdz_t_slice_const_lat_const_depth));
% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lat_idx_varid, lat_idx);
netcdf.putVar(ncid, lat_varid, lat_rho_vec(lat_idx));
netcdf.putVar(ncid, lon_varid, lon_rho_vec(:));
netcdf.putVar(ncid, w_varid, w_t_slice_const_lat_const_depth);
netcdf.putVar(ncid, temp_varid, temp_t_slice_const_lat_const_depth);
netcdf.putVar(ncid, salt_varid, salt_t_slice_const_lat_const_depth);
netcdf.putVar(ncid, dTdz_varid, dTdz_t_slice_const_lat_const_depth);
% close netcdf file
netcdf.close(ncid);
%----------------------------------------
% save data into a netcdf file
disp('Done saving w_temp_salt_slice_const_lat_cost_depth!')
%--------------------------------------------------------
