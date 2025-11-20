%

% save timeseries of box avg vorticity for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
%%Loop to calculate averages if need (modify as needed)
indxRange = 0:3877; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
time_arr = zeros(Nt, 1);
f0 = 0.909*1e-4; % Coriolis freq ; 
%% plot options
vlevel_skip = 25;
legendCell = [];
% vlevel_arr = ([100, 52, 39, 15, 1]);% NumLayers: -vlevel_skip: 1
vlevel = -500; % vlevel < 0 => depth;
Nz = 1; % length(vlevel_arr);
vort_t = zeros(Nt, Ny, Nx);
box_avg_vort_t = zeros(Nt, 1);
box="AtlantisII_1km" % average around AtlantisII
% box = "domain_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min

% read speed_t_hslice
filename = strcat('sm3_2019_2020_vort_t_hslice_vlevel_', string(vlevel) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
vort_t_hslice = ncread(filename, 'vort_t_hslice');
t_arr = ncread(filename, 't_arr');
%---------------------------------------
% Time loop
%
box_avg_vort_t = zeros(Nt, 1);

for nt=1:Nt
    vort_hslice_box = squeeze(vort_t_hslice(nt, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    
    vort_hslice_box = vort_hslice_box./f0;

    box_avg_vort_t(nt) = mean(vort_hslice_box, 'all', 'omitnan');
    sprintf('box_avg_vort_t = %.7e', box_avg_vort_t(nt) )

end
%-------------------------------------------
%%
% save data into a netcdf file

if(strcmp(box, 'AtlantisII_1km')) % covers only Atlantis II
        vort_t_filename = strcat('sm3_2019_2020_box_avg_vort_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');                                                                                               % ncid = netcdf.create(filename,'CLOBBER');
elseif(strcmp(box, 'NESM_1km')) % covers all three seamounts
        vort_t_filename = strcat('sm3_2019_2020_bigbox_avg_vort_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
        % ncid = netcdf.create(filename,'CLOBBER');
elseif(strcmp(box, 'domain_1km')) % covers all three seamounts
        vort_t_filename = strcat('sm3_2019_2020_domain_avg_vort_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
end
ncid = netcdf.create(vort_t_filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(vort_t_filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');


vort_varid = netcdf.defVar(ncid, 'box_avg_vort_t', 'double', [t_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'timeseries of (dimless) box avg vorticity (zeta/f) at a given vlevel');
netcdf.putAtt(ncid, vort_varid, 'units', '--');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(box_avg_vort_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, vort_varid, box_avg_vort_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating box_avg_vort_t for selected depth layers!')

%%

   
