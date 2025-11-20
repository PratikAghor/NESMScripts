%
% save time series of vertical slices of ertel potential vorticity for a selected lat values
% 
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
% indxRange = 952:3877; % entire year
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
indxRange=3024:3152; % Jan 15 - Jan 30, 2020
% indxRange=3264:3392; % Feb 14 - Feb 29, 2020

init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);

legendCell = [];
lat_idx_arr = ([90, 112, 135]);
lat_arr = zeros(length(lat_idx_arr), 1);
for i = 1:length(lat_arr_idx)
	lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
lat_arr
lonsec = lon_rho_vec;

pv_3d_t = zeros(Nt, NumLayers, Ny, Nx);
rho_pot_3d_t = zeros(Nt, NumLayers, Ny, Nx);

sprintf('size(pv_3d_t) = %d', size(pv_3d_t))
%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_avg.%05d.nc'], ' file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    zeta = pagetranspose(ncread(hisfile, 'zeta'));
    time = ncread(fname, 'time');
    time_arr(nt) = time;
    % method 1
    epv = ertel_aghor(hisfile, gridfile, 'rho', tindex); % get ertel pv, 3d mat
    % size(epv)
    epv = epv(2:end, :, :);
    % epv = shiftdim(epv, 1);
    % max((epv(40, :, :)))
    % size(epv)
    pv_3d_t(nt, :, :, :) = -g.*epv;
    
    % save rho_pot as well
    % get 3d temp, salt in (Nz, Ny, Nx) form
    temp = pagetranspose(ncread(hisfile, 'temp'));
    salt = pagetranspose(ncread(hisfile, 'salt'));

    temp = shiftdim(temp, 2);
    salt = shiftdim(salt, 2);

    rhotmp = gsw_sigma2(salt, temp); % get potential density
    size(rhotmp)

    rho_pot_3d_t(nt, :, :, :)  = rhotmp; % get potential density 

end

%%
% save data into a netcdf file

vort_t_filename = strcat('nesm_2019_2020_pv_3d_t', '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
ncid = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
depth_len = netcdf.defDim(ncid, 'Nz', NumLayers);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

% time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
% netcdf.putAtt(ncid, time_varid, 'description', 'time array');
% netcdf.putAtt(ncid, time_varid, 'units', 's');

% vlevel_varid = netcdf.defVar(ncid, 'vlevel', 'double', [depth_len]);
% netcdf.putAtt(ncid, vlevel_varid, 'description', 'array of vertical layers from the original data');
% netcdf.putAtt(ncid, vlevel_varid, 'units', '--');
% netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(depth_len));

vort_varid = netcdf.defVar(ncid, 'pv_3d_t', 'double', [t_len depth_len y_len x_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of 3d ertel potential vorticity');
netcdf.putAtt(ncid, vort_varid, 'units', 's^-3');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(pv_3d_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
% netcdf.putVar(ncid, vlevel_varid, vlevel_arr);
netcdf.putVar(ncid, vort_varid, pv_3d_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating pv_3d_t!')

%%
%------------------------------------------------------------------
% save data into a netcdf file

vort_t_filename = strcat('nesm_2019_2020_rho_pot_3d_t', '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
ncid = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
depth_len = netcdf.defDim(ncid, 'Nz', NumLayers);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

rho_varid = netcdf.defVar(ncid, 'rho_pot_3d_t', 'double', [t_len depth_len y_len x_len]);
netcdf.putAtt(ncid, rho_varid, 'description', 'time series of 3d potential density at 2000 dbar from TEOS10');
netcdf.putAtt(ncid, rho_varid, 'units', 'kgm^-3');
netcdf.putAtt(ncid, rho_varid, 'array dimensions', size(rho_pot_3d_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
% netcdf.putVar(ncid, vlevel_varid, vlevel_arr);
netcdf.putVar(ncid, rho_varid, rho_pot_3d_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating rho_pot_3d_t!')
%----------------------------------------------------------------
