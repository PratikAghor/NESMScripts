% assume annual mean of ub, vb, wb has been saved.
% calculate up = u - ub, vp = v - vb, etc.
% calculate eke = 0.5*(up^2 + vp^2)
% convert the eke from sigma levels to zlevs
% take water column avg from z = -1800 to z = -5000
% save the timeseries
%------------------------------------------------
start_paths
[~, ~, mask] = read_latlonmask(gridfile, 'r');
%------------------------------------------------
indxRange=0:3877; % entire year
[~, Nt] = size(indxRange);

% annual_mean_indxRange = 952:3877;
% [~, mean_Nt] = size(annual_mean_indxRange);
% annual_mean_uvw_file = strcat('nesm_2019_2020_uvw_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

% ub = ncread(annual_mean_uvw_file, 'ub');
% vb = ncread(annual_mean_uvw_file, 'vb');
% wb = ncread(annual_mean_uvw_file, 'wb');
% zetab = ncread(annual_mean_uvw_file, 'zetab');
%------------------------------------------------
vlevel_arr = ([-3000]); % vlevel < 0 => depth;
Nz = length(vlevel_arr);

box="AtlantisII_1km" % average around AtlantisII
% box="NESM_1km"
% box = "domain_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%------------------------------------------------
% read speed_t_hslice
filename = strcat('nesm_2019_2020_eke_t_hslice_vlevel_', string(vlevel_arr(1)) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
eke_t_hslice = ncread(filename, 'eke_t_hslice');
t_arr = ncread(filename, 't_arr');
%---------------------------------------
% Time loop
%
box_avg_eke_t = zeros(Nt, 1);

for nt=1:Nt
    eke_hslice_box = squeeze(eke_t_hslice(nt, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    
    % disp('flag')
    if(nt == 715)
	    eke_hslice_box
    end

    box_avg_eke_t(nt) = mean(eke_hslice_box, 'all', 'omitnan');
    sprintf('box_avg_eke_t = %.7e', box_avg_eke_t(nt) )

end
% size(eke_t_box_avg)
%----------------------------------------------------------------------------
% save data in netcdf
%--------------------------------------------------------
if(strcmp(box, 'AtlantisII_1km')) % covers only Atlantis II
        filename = strcat('nesm_2019_2020_box_avg_eke_t_hslice_vlevel_', string(vlevel_arr(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
	% ncid = netcdf.create(filename,'CLOBBER');
elseif(strcmp(box, 'NESM_1km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_bigbox_avg_eke_t_hslice_vlevel_', string(vlevel_arr(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
	% ncid = netcdf.create(filename,'CLOBBER');
elseif(strcmp(box, 'domain_1km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_domain_avg_eke_t_hslice_vlevel_', string(vlevel_arr(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
end

ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');
netcdf.putAtt(ncid, time_varid, 'array dimensions', size(t_arr));

vort_varid = netcdf.defVar(ncid, 'box_avg_eke_t', 'double', [t_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of <EKE>');
netcdf.putAtt(ncid, vort_varid, 'units', 'm^2/s^2');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(box_avg_eke_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, t_arr); % start from 0
netcdf.putVar(ncid, vort_varid, box_avg_eke_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating box_avg_eke_t_hslice for selected depth layers!')
%--------------------------------------------------------
