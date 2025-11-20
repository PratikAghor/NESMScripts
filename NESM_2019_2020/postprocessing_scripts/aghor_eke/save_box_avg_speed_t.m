%
% save box averaged speed t, assume speed_t_hslice has been saved
% and time-stamps
%
%-------------------------------------------------------
start_paths
type='r'; % see get_hslice
% tindex=1; 
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange = 0:3877; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);

%% plot options
vlevel_skip = 25;
legendCell = [];
% vlevel_arr = ([100, 52, 39, 15, 1]);% NumLayers: -vlevel_skip: 1
vlevel_arr = ([-3000]); % vlevel < 0 => depth;
Nz = length(vlevel_arr);
w_t = zeros(Nt, Ny, Nx, length(vlevel_arr));
box_avg_speed_t = zeros(Nt, 1);
box="AtlantisII_1km" % average around AtlantisII
% box="NESM_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%---------------------------------------
% read speed_t_hslice
filename = strcat('nesm_2019_2020_speed_t_hslice_vlevel_', string(vlevel_arr(1)) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
speed_t_hslice = ncread(filename, 'u_t');
t_arr = ncread(filename, 't_arr');
%---------------------------------------

%%
for nt = 1:Nt    
    box_speed_hslice = squeeze(speed_t_hslice(nt, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    box_avg_speed_t(nt) = squeeze(nanmean(nanmean(box_speed_hslice, 1), 2));    
end

%%
% save data into a netcdf file
if(strcmp(box, 'AtlantisII_1km')) % covers only Atlantis II
        filename = strcat('nesm_2019_2020_box_avg_speed_t_hslice_vlevel_', string(vlevel_arr(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');ncid = netcdf.create(filename,'CLOBBER');
elseif(strcmp(box, 'NESM_1km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_bigbox_avg_speed_t_hslice_vlevel_', string(vlevel_arr(1)), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');ncid = netcdf.create(filename,'CLOBBER');
end

% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
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

vlevel_varid = netcdf.defVar(ncid, 'vlevel', 'double', [depth_len]);
netcdf.putAtt(ncid, vlevel_varid, 'description', 'array of vertical layers from the original data');
netcdf.putAtt(ncid, vlevel_varid, 'units', '--');
netcdf.putAtt(ncid, vlevel_varid, 'array dimensions', size(depth_len));

vort_varid = netcdf.defVar(ncid, 'box_avg_speed_t', 'double', [t_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of <speed>');
netcdf.putAtt(ncid, vort_varid, 'units', 'm/s');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(box_avg_speed_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, vlevel_varid, vlevel_arr);
netcdf.putVar(ncid, vort_varid, box_avg_speed_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating box_avg_speed_t_hslice for selected depth layers!')

%%

   
