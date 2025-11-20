%
% save horizontal slices of w for selected vertical levles (vlevels)
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

vlevel_skip = 25;
legendCell = [];

vlevel = -4000; % vlevel < 0 => depth;
Nz = 1; % length(vlevel_arr);
pv_t = zeros(Nt, Ny, Nx);
box_avg_pv_t = zeros(Nt, 1);
box="AtlantisII_5km" % average around AtlantisII
% box="NESM_1km"
% box="domain_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%--------------------------------------------------
% read speed_t_hslice
filename = strcat('nesm_2019_2020_5km_epv_t_hslice_vlevel_', string(vlevel) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
pv_t_hslice = ncread(filename, 'epv_t_hslice');
t_arr = ncread(filename, 't_arr');
%---------------------------------------
%%
for nt = 1:Nt
    % sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
    % fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
    % hisfile = fname;
    % time = ncread(fname, 'time');
    % time_arr(nt) = time;
    % for k = 1:Nz % index for layer numbers (since we don't save all layers)
    %    vlevel = vlevel_arr(k); % actual layer number
    %    w_t(nt, :, :, k)=get_hslice(hisfile,gridfile,'w',tindex,vlevel,type); 
    % end
    box_pv_t = squeeze(pv_t_hslice(nt, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max));
    box_avg_pv_t(nt) =  mean(box_pv_t, 'all', 'omitnan');   
    sprintf('box_avg_pv_t = %.7e', box_avg_pv_t(nt) )	
end

%%
% save data into a netcdf file
if(strcmp(box, 'AtlantisII_5km')) % covers only Atlantis II
        filename = strcat('nesm_2019_2020_5km_box_avg_pv_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
elseif(strcmp(box, 'NESM_5km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_5km_bigbox_avg_pv_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
elseif(strcmp(box, 'domain_5km')) % covers all three seamounts
        filename = strcat('nesm_2019_2020_5km_domain_avg_pv_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
end

ncid = netcdf.create(filename,'CLOBBER');
% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% depth_len = netcdf.defDim(ncid, 'Nz', Nz);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');

vort_varid = netcdf.defVar(ncid, 'box_avg_pv_t', 'double', [t_len]);
netcdf.putAtt(ncid, vort_varid, 'description', 'time series of w hslice');
netcdf.putAtt(ncid, vort_varid, 'units', 's^-1');
netcdf.putAtt(ncid, vort_varid, 'array dimensions', size(box_avg_pv_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, time_arr); % start from 0
% netcdf.putVar(ncid, lat_varid, lat_rho);
% netcdf.putVar(ncid, lon_varid, lon_rho);
netcdf.putVar(ncid, vort_varid, box_avg_pv_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating box_avg_pv_t_hslice for selected depth layers!')

%%

   
