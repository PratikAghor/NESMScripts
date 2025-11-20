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
indxRange=3016:3072; % Jan 14 - Jan 20, 2020
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
%------------------------------------------------------
% for vslice at const_lat
lat_idx_arr = ([112]);

lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
for i = 1:length(lat_idx_arr)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lonsec = lon_rho_vec;
X = (repmat(lon_rho_vec, N, 1));
%------------------------------------------------------
% for vslice at const lon
lon_idx_arr = ([172]);
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));
end
latsec = lat_rho_vec;
%---------------------------------------------------

%---------------------------------------------------
vort_t_vslice_const_lon = zeros(Nt, N, Ny, length(lon_idx_arr));
t_arr = zeros(Nt, 1);

for nt=1:Nt
	sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
	fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
	hisfile = fname;
	
	t_arr(nt) = ncread(hisfile, 'time');

    	vort_t_slices = zeros(NumLayers, Ny, Nx);

    	% Calculate vorticity for each layer
    	for vlevel = 1:NumLayers
        	[~, ~, ~, vort_t_slices(vlevel, :, :)]=get_vort(hisfile,gridfile,tindex,vlevel,coef);
    	end

    	% for k = 1:length(lat_arr) % index for layer numbers
        %	f0 = 2*Omega*sin(lat_arr(k)); % Coriolis freq.
        %	lat_val = lat_arr(k); % actual const latitude
        %	lat_idx = find(abs(lat_rho_vec-lat_arr(k))<1e-3); % find idx of lat in lat_rho_vec
        %	% Normalize with coriolis freq.
        %	vort_t_vslice(nt, k, :, :) = squeeze(vort_t_slices(lat_idx, :, :))./f0;
    	% end

    	for j =1:length(lon_arr)
        	lon_idx = lon_idx_arr(j);
        	lon_val = lon_arr(j); % actual const longitude
        	vort_t_vslice_const_lon(nt, :, :, j) = squeeze(vort_t_slices(:, :, lon_idx));
    	end

end
%--------------------------------------
% save data into a netcdf file

filename = strcat('nesm_2019_2020_vort_t_vslice_const_lon_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', length(lon_idx_arr));
y_len = netcdf.defDim(ncid, 'Ny', Ny);
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

lon_idx_varid = netcdf.defVar(ncid, 'lon_idx_arr', 'double', [x_len]);
netcdf.putAtt(ncid, lon_idx_varid, 'description', 'indices in lon_rho_vec for const lon values where vertical slice is taken');
netcdf.putAtt(ncid, lon_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lon_idx_varid, 'array dimensions', size(lon_idx_arr));

lon_varid = netcdf.defVar(ncid, 'lon_arr', 'double', [x_len]);
netcdf.putAtt(ncid, lon_varid, 'description', 'const lon values where vertical slice is taken');
netcdf.putAtt(ncid, lon_varid, 'units', '--');
netcdf.putAtt(ncid, lon_varid, 'array dimensions', size(lon_arr));

w_varid = netcdf.defVar(ncid, 'vort_t_vslice_const_lon', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'timeseries of vort_t at const_lon');
netcdf.putAtt(ncid, w_varid, 'units', 's^{-1}');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(vort_t_vslice_const_lon));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lon_idx_varid, lon_idx_arr);
netcdf.putVar(ncid, lon_varid, lon_arr);
netcdf.putVar(ncid, w_varid, vort_t_vslice_const_lon);
% close netcdf file
netcdf.close(ncid);
%----------------------------------------
% save data into a netcdf file
disp('Done saving vort_t_vslice_const_lon!')
%--------------------------------------------------------
