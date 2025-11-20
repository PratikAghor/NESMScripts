% assume rho_w_annual_mean has been saved
% Calculate VEBF = \bar{w'b'} with respect to the annual mean
% Pratik Aghor
%%
%---------------------------------------------------
clc
clear all
close all
start_paths
%---------------------------------------------------
% indxRange = 469:3382; % entire year
% indxRange = 2080:2144; % Sept 2019
indxRange = 3016:3072; % 14 Jan-20 Jan 2020
% indxRange = 3080:3144; Jan 2020

nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
%---------------------------------------------------
annual_mean_indxRange = 952:3877;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_file = strcat('nesm_2019_2020_5km_rho_w_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

wb = ncread(annual_mean_file, 'wb');
rhob = ncread(annual_mean_file, 'rhob');

% size(wb)
% size(rhob)
%---------------------------------------------------
% wprime = zeros(N, M, L);
% bprime = zeros(N, M, L);
VEBF3D_avg = zeros(N, M, L);

for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;

    rhotmp = zeros(N, M, L);
    wtmp = zeros(N, M, L);

    % construct 3d rho mat from get_rho.m
    for k=1:N
        vlevel = k;
        % tmp = get_rho(hisfile, gridfile, tindex, vlevel, coef);
        % size(tmp)
        [~, ~, ~, rhotmp(k, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    end
    wtmp = pagetranspose(ncread(hisfile, 'w'));
    wtmp = shiftdim(wtmp, 2);
    size(wtmp)
    size(rhotmp)
 
    wprime = wtmp - wb;
    bprime =(-g/rho0).*(rhotmp - rhob);
   
    wpbp = wprime.*bprime;
 
    VEBF3D_avg = VEBF3D_avg + (wpbp);
end
VEBF3D_avg = VEBF3D_avg./Nt;
%---------------------------------------------------
%
% Save
%
% save data into a netcdf file
filename = strcat('nesm_2019_2020_5km_VEBF_3d_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
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

vebf_varid = netcdf.defVar(ncid, 'VEBF3D', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, vebf_varid, 'description', '3d VEBF (wprime*bprime) time avg (Nz, Ny, Nx)');
netcdf.putAtt(ncid, vebf_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vebf_varid, 'array dimensions', size(VEBF3D_avg));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, vebf_varid, VEBF3D_avg);
% close netcdf file
netcdf.close(ncid);
disp('Done creating VEBF3D_avg!')
%--------------------------------------------------------
