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

annual_mean_indxRange = 952:3877;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_uvw_file = strcat('nesm_2019_2020_uvw_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

ub = ncread(annual_mean_uvw_file, 'ub');
vb = ncread(annual_mean_uvw_file, 'vb');
% wb = ncread(annual_mean_uvw_file, 'wb');
% zetab = ncread(annual_mean_uvw_file, 'zetab');
%------------------------------------------------
vlevel_arr = ([-3000]); % vlevel < 0 => depth;
Nz = length(vlevel_arr);

box="AtlantisII_1km" % average around AtlantisII
% box="NESM_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%------------------------------------------------
% Time loop
%
t_arr = zeros(Nt, 1);
box_avg_eke_t = zeros(Nt, 1);
eke_t_hslice = zeros(Nt, Ny, Nx);
for nt=1:Nt
   
    sprintf('nt = %05d', indxRange(1, nt))

    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    
    u = pagetranspose(ncread(hisfile, 'u'));
    v = pagetranspose(ncread(hisfile, 'v'));
    u = shiftdim(u, 2);
    v = shiftdim(v, 2);

    t_arr(nt) = ncread(hisfile, 'time');

    % zeta = pagetranspose(ncread(hisfile, 'zeta'));
    up = u - ub; % u' = u - \bar{u}
    vp = v - vb; % v' = v - \bar{v}	
   
    
    up = u2rho_3d(up);
    vp = v2rho_3d(vp);
    size(up)
    size(vp)
    % up_3d_zlevs = zeros(Nz, Ny, Nx);
    % vp_3d_zlevs = zeros(Nz, Ny, Nx);

    % get hslices of up, vp at zlevs
    
	vlevel = vlevel_arr(1);
	up_hslice = vinterp(up, zr, vlevel);
        vp_hslice = vinterp(vp, zr, vlevel);

	if(nt == 10)
		up_hslice
	end
        eke_t_hslice(nt, :, :) = 0.5.*(up_hslice.^2 + vp_hslice.^2); % assignment requires large memory, apparently
        
        
     
    % eke_hslice_box = squeeze(eke_t_hslice(nt, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    
    disp('flag')
    % size(eke_hslice_box)
    % box_avg_eke_t(nt) = mean(eke_hslice_box, 'all', 'omitnan');
end
% size(eke_t_box_avg)
%----------------------------------------------------------------------------
% save
% save data into a netcdf file
filename = strcat('nesm_2019_2020_eke_t_hslice_vlevel_', string(vlevel_arr(1)) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');

ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');
netcdf.putAtt(ncid, time_varid, 'array dimensions', size(t_arr));

eke_varid = netcdf.defVar(ncid, 'eke_t_hslice', 'double', [t_len y_len x_len]);
netcdf.putAtt(ncid, eke_varid, 'description', 'EKE hslice at const vlevel (Nt)');
netcdf.putAtt(ncid, eke_varid, 'units', 'm^2s^-2');
netcdf.putAtt(ncid, eke_varid, 'array dimensions', size(eke_t_hslice));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, t_arr); % start from 0
netcdf.putVar(ncid, eke_varid, eke_t_hslice);
% close netcdf file
netcdf.close(ncid);
disp('Done creating eke_t_hslice at const vlevel!')
%--------------------------------------------------------
