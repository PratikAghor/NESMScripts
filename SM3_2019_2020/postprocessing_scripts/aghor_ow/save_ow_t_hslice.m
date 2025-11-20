%
% save horizontal slices of vorticity for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange=952:3877; % entire year
% indxRange=2032:2160; % Sep 15 - Sep 30, 2019
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
% indxRange=3024:3152; % Jan 15 - Jan 30, 2020

init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
t_arr = zeros(Nt, 1);

% vlevel_arr = ([100, 52, 39, 15, 1]);% NumLayers: -vlevel_skip: 1
vlevel = -4000; % vlevel < 0 => depth;
Nz = 1; % length(vlevel_arr);
ow_t = zeros(Nt, Ny, Nx);
%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'SM3_2019_2020_avg.%05d.nc'], ' file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'SM3_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    t_arr(nt) = time;
    [lat,lon,mask, ow_t(nt, :, :)]=get_okubo(hisfile,gridfile,tindex,vlevel,coef);   
end

%%
% save data into a netcdf file

ow_t_filename = strcat('sm3_2019_2020_ow_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
ncid = netcdf.create(ow_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(ow_t_filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

vort_varid = netcdf.defVar(ncid, 'ow_t_hslice', 'double', [t_len y_len x_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'Okubo-Weiss parameter hslice at different vlelevs');
netcdf.putAtt(ncid, vort_varid, 'units', 's^-2');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(ow_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, t_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, vort_varid, ow_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating ow_t for selected depth layers!')

%%

   
