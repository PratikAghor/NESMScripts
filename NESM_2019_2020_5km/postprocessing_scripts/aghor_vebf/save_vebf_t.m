% assume rho_w_annual_mean has been saved
% Calculate VEBF = \bar{w'b'} with respect to the annual mean
% save timeseries of a horizontal avg at a const depth
% Pratik Aghor
%%
%---------------------------------------------------
clc
clear all
close all
start_paths;
%---------------------------------------------------
% box limits for horizontal averages
box="AtlantisII_5km" % average around AtlantisII
% box = "domain_5km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%---------------------------------------------------
indxRange = 952:3877; % entire year
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
vlevel = -4000;
%---------------------------------------------------
annual_mean_indxRange = 952:3877;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_file = strcat('nesm_2019_2020_5km_rho_w_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

wb = ncread(annual_mean_file, 'wb');
rhob = ncread(annual_mean_file, 'rhob');

% size(wb)
% size(rhob)
%---------------------------------------------------
vebf_t_box_avg = zeros(Nt, 1);

for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    %----------------------------------------------------
    % method 1: calculate VEBF 3d and then take an hslice, slow
    % rhotmp = zeros(N, M, L);
    % wtmp = zeros(N, M, L);

    % construct 3d rho mat from get_rho.m
    % for k=1:N
    %    vlevel_tmp = k;
    %    % tmp = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    %    % size(tmp)
    %    [~, ~, ~, rhotmp(k, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel_tmp, coef);
    % end
    % wtmp = pagetranspose(ncread(hisfile, 'w'));
    % wtmp = shiftdim(wtmp, 2);
    % size(wtmp)
    % size(rhotmp)

    % wprime = wtmp - wb;
    % bprime =(-g/rho0).*(rhotmp - rhob);

    % wpbp = wprime.*bprime;

    % VEBF3D = (wpbp);
    
    % % updated from save_KmKe_2d.m, take an hslice at a given vlevel
    % zeta = pagetranspose(ncread(hisfile, 'zeta'));

    % zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
    % zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);

    % VEBF2D = vinterp(VEBF3D,zr,vlevel);
    %----------------------------------------------------
    % method 2: calculate rho_hslice w_hslice and VEBF2D directly, fast
    wtmp = pagetranspose(ncread(hisfile, 'w'));
    wtmp = shiftdim(wtmp, 2);
    size(wtmp)

    w_hslice = vinterp(wtmp, zr, vlevel);
    wb_hslice = vinterp(wb, zr, vlevel);
    wprime = w_hslice - wb_hslice;

    [~, ~, ~, rho_hslice] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    rhob_hslice = vinterp(rhob, zr, vlevel);
    rhoprime = rho_hslice - rhob_hslice;

    bprime = (-g/rho0).*rhoprime;

    VEBF2D = wprime.*bprime;
    %----------------------------------------------------
    % take a horizontal avg around the middle seamount
    VEBF_t_box = VEBF2D(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);

    VEBF_t_box_avg(nt) = mean(VEBF_t_box, 'all', 'omitnan');
    
    VEBF_t_box_avg(nt)    
end
%---------------------------------------------------
%---------------------------------------------------
%
% Save
%
% save data into a netcdf file
if(strcmp(box, 'AtlantisII_5km')) % covers only Atlantis II
        filename = strcat('nesm_2019_2020_5km_box_avg_vebf_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', string(indxRange(Nt)),'.nc');

elseif(strcmp(box, 'domain_5km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_5km_box_avg_vebf_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', string(indxRange(Nt)),'.nc');
end
ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
z_len = netcdf.defDim(ncid, 'Nz', N);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

vebf_varid = netcdf.defVar(ncid, 'VEBF_t', 'double', [t_len]);
netcdf.putAtt(ncid, vebf_varid, 'description', '3d VEBF (wprime*bprime) time series at a const depth (Nt)');
netcdf.putAtt(ncid, vebf_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vebf_varid, 'array dimensions', size(VEBF_t_box_avg));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vebf_varid, VEBF_t_box_avg);
% close netcdf file
netcdf.close(ncid);
disp('Done creating VEBF_t!')
%--------------------------------------------------------
