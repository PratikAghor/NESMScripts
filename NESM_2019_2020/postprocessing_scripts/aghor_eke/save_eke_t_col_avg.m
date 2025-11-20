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
indxRange=469:3382; % entire year
[~, Nt] = size(indxRange);

annual_mean_indxRange = 469:3382;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_uvw_file = strcat('nesm_2019_2020_uvw_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

ub = ncread(annual_mean_uvw_file, 'ub');
vb = ncread(annual_mean_uvw_file, 'vb');
% wb = ncread(annual_mean_uvw_file, 'wb');
% zetab = ncread(annual_mean_uvw_file, 'zetab');
%------------------------------------------------
% save 50 zlevs
Nz = 50; 
zlevs = linspace(-1800, -5000, N);
vlevel_arr = zlevs;
%------------------------------------------------
% Time loop
%
eke_t_col_avg = zeros(Nt, 1);
eke_3d_zlevs = zeros(Nt, Nz, Ny, Nx);
size(eke_3d_zlevs)
for nt=1:Nt
    
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    
    u = pagetranspose(ncread(hisfile, 'u'));
    v = pagetranspose(ncread(hisfile, 'v'));
    u = shiftdim(u, 2);
    v = shiftdim(v, 2);
    % zeta = pagetranspose(ncread(hisfile, 'zeta'));
    up = u - ub; % u' = u - \bar{u}
    vp = v - vb; % v' = v - \bar{v}	
   
    
    up = u2rho_3d(up);
    vp = v2rho_3d(vp);
    size(up)
    size(vp)
    up_3d_zlevs = zeros(Nz, Ny, Nx);
    vp_3d_zlevs = zeros(Nz, Ny, Nx);

    % get hslices of up, vp at zlevs
    for k = 1:Nz
	vlevel = vlevel_arr(k);
	up_3d_zlevs(k, :, :) = vinterp(up, zr, vlevel);
        vp_3d_zlevs(k, :, :) = vinterp(vp, zr, vlevel);

        up_hslice = squeeze(up_3d_zlevs(k, :, :));
	vp_hslice = squeeze(vp_3d_zlevs(k, :, :));
	if(k==1)
		size(up_hslice)
		size(vp_hslice)
	end	
	eke_3d_zlevs(nt, k, :, :) = 0.5.*(up_hslice.^2 + vp_hslice.^2); % assignment requires large memory, apparently
        
        % eke_3d_zlevs(nt, k, :, :) = eke; 
    end
     
    eke_3d_zlevs_box = squeeze(eke_3d_zlevs(nt, :, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    
    disp('flag')
    size(eke_3d_zlevs_box)
    eke_t_col_avg(nt) = mean(eke_3d_zlevs_box, 'all', 'omitnan');
end
size(eke_t_col_avg)
%----------------------------------------------------------------------------
% save
% save data into a netcdf file
filename = strcat('nesm_2019_2020_eke_t_col_avg_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

eke_varid = netcdf.defVar(ncid, 'eke_col_avg', 'double', [t_len]);
netcdf.putAtt(ncid, eke_varid, 'description', 'EKE column average within box near Atlantis seamount (Nt)');
netcdf.putAtt(ncid, eke_varid, 'units', 'm^2s^-2');
netcdf.putAtt(ncid, eke_varid, 'array dimensions', size(eke_t_col_avg));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, eke_varid, eke_t_col_avg);
% close netcdf file
netcdf.close(ncid);
disp('Done creating EKE2D_col_avg!')
%--------------------------------------------------------
