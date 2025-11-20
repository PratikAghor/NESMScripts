%
% save horizontal slices of w for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
type='r'; % see get_hslice
% tindex=1; 
%-------------------------------------------------------
indxRange = 952:3877; % entire year

% z lims, z1 < z2
z1 = -4000;
z2 = -1000; % -500;

% Aghor's script
q_thresh = 1e-9; % threshold for low PV [s^-3]

Nt = numel(indxRange);

t_arr = zeros(Nt, 1);
fq = zeros(Nt,1); % fraction of volume with |q| < q_thresh;
z_box_avg_pv_t = zeros(Nt,1); % V_q
t_arr = zeros(Nt,1);
IS_V = zeros(Nt, 1); % instability strength

box = "AtlantisII_1km";
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);


for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SM3_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SM3_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    t_arr(nt) = time;
    epv_3d = ertel_aghor(hisfile,gridfile,'rho',tindex);
    epv_3d = -g.*epv_3d;

    size(epv_3d)

    % Interpolate epv to rho levels
    epv_rho = 0.5 * (epv_3d(1:N, :, :) + epv_3d(2:N+1, :, :));  % (N x Ny x Nx)

    size(epv_rho)

    % Identify indices corresponding to 1000–4000 m
    z_mask = (zr <= z2) & (zr >= z1);  % zr is (N x Ny x Nx)

    % Apply mask to horizontal box and depth range
    box_mask = false(size(epv_rho));
    box_mask(:, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max) = true;
    full_mask = box_mask & z_mask;

    % Extract EPV values in the full 3D prism
    box_epv_vals = epv_rho(full_mask);

    % Compute fraction of |q| < threshold
    validMask = ~isnan(box_epv_vals);
    fq(nt) = sum(abs(box_epv_vals(validMask)) < q_thresh) / sum(validMask);

    % vol avg of q in the horizontal box around Atlantis II bounedd b/w z1 and z2    
    z_box_avg_pv_t(nt) = mean(box_epv_vals(validMask), 'all', 'omitnan');


    IS_V(nt) = fq(nt) .* abs(z_box_avg_pv_t(nt));  % instantaneous contribution
 
    fprintf('fq = %.3g, z_box_avg_pv_t = %.3g\n', fq(nt), z_box_avg_pv_t(nt))

    sprintf("IS_V(nt) = %.2f", IS_V(nt)) 
end

%%
%----------------------------------------------------------------------------
% save
% save data into a netcdf file
filename = sprintf('sm3_2019_2020_IS_V_t_z1_%d_z2_%d_nt_%d_%d.nc', ...
    z1, z2, indxRange(1), indxRange(Nt));

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

% fq
fq_varid = netcdf.defVar(ncid, 'fq', 'double', [t_len]);
netcdf.putAtt(ncid, fq_varid, 'description', 'fraction of low PV (|q| < q_thresh) points');
netcdf.putAtt(ncid, fq_varid, 'units', 'dimensionless');
netcdf.putAtt(ncid, fq_varid, 'array dimensions', size(fq));

% z_box_avg_pv_t
zavg_varid = netcdf.defVar(ncid, 'z_box_avg_pv_t', 'double', [t_len]);
netcdf.putAtt(ncid, zavg_varid, 'description', 'box-averaged volume-mean PV');
netcdf.putAtt(ncid, zavg_varid, 'units', 's^-3');
netcdf.putAtt(ncid, zavg_varid, 'array dimensions', size(z_box_avg_pv_t));

% IS_V
isv_varid = netcdf.defVar(ncid, 'IS_V', 'double', [t_len]);
netcdf.putAtt(ncid, isv_varid, 'description', 'Instantaneous IS_V = fq * |<q>_V|');
netcdf.putAtt(ncid, isv_varid, 'units', 's^-3');
netcdf.putAtt(ncid, isv_varid, 'array dimensions', size(IS_V));

%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, t_arr); % start from 0
netcdf.putVar(ncid, fq_varid, fq);
netcdf.putVar(ncid, zavg_varid, z_box_avg_pv_t);
netcdf.putVar(ncid, isv_varid, IS_V);
% close netcdf file
netcdf.close(ncid);
disp('Done creating IS_V time series!')
%--------------------------------------------------------
%%

   
