% save_xyz_avg_kv_t.m
%--------------------------------------------------------------------------
% Save vertically (z=-1800 to zbottom) and box averaged Kv = <Kv> versus t
% from -1800 m to bottom, plus horizontal box averaging
% Author: Pratik Aghor
%--------------------------------------------------------------------------

clear all; close all; clc;
start_paths;  % Your setup script that defines paths, zw, etc.

indxRange = 952:3877; % entire year
Nt = length(indxRange);

vlevel_top = -1000;  % top depth of averaging layer
% Bottom is model bottom (assumed at zw(end))

Kmax = 5e-3; % max kv value

domain_avg_kv_t = zeros(Nt,1);
box_avg_kv_t = zeros(Nt,1);
t_arr = zeros(Nt,1);

box = "AtlantisII_1km"; % box for horizontal averaging
% box = "domain_1km"; % bigger domain box

[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);

% Find vertical indices corresponding to layer between vlevel_top and bottom
% bottom_depth = -5000; % zw(1); % zw goes from bottom height to top
% idx_layer = find(zw <= vlevel_top & zw >= bottom_depth); 

for nt = 1:Nt
    fname = sprintf([dirHR, 'SM3_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);

    t_arr(nt) = ncread(fname, 'time');

    AKt = ncread(fname, 'AKt'); % size [Nx, Ny, Nz+1]
    AKt = shiftdim(AKt, 2); % now it is (Nz+1) x Ny x Nx, like I prefer, zw also has (Nz+1) x Ny x Nx size

     AKt(AKt > Kmax) = Kmax; % clip to max value
    if(nt == 1)
		size(AKt)
    end
   
    % Create a mask where z_w < vlevel_top (i.e., deeper than 1800 m)
    mask_layer = zw <= vlevel_top;           % [Nz+1 x Ny x Nx]

    % Apply mask to AKt
    AKt_layer = AKt;
    AKt_layer(~mask_layer) = NaN;

    % Domain average over full area
    domain_avg_kv_t(nt) = mean(AKt_layer(:), 'omitnan');

    % Box average
    AKt_box = AKt_layer(:, lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);
    box_avg_kv_t(nt) = mean(AKt_box(:), 'omitnan');

    
    fprintf('t=%d: domain avg Kv = %.4e, box avg Kv = %.4e\n', nt, domain_avg_kv_t(nt), box_avg_kv_t(nt));
end

%% Save to NetCDF

% if strcmp(box, 'AtlantisII_1km')
filename = sprintf('sm3_2019_2020_model_xyz_avg_kv_t_%d_to_bottom_nt_%d_%d.nc', round(vlevel_top), indxRange(1), indxRange(end));

% Define dimensions

ncid = netcdf.create(filename, 'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);

% Define variables and attributes
t_varid = netcdf.defVar(ncid, 't_arr', 'double', t_len);
netcdf.putAtt(ncid, t_varid, 'description', 'time array measured from 01/01/2018');
netcdf.putAtt(ncid, t_varid, 'units', 's');

domain_avg_kv_varid = netcdf.defVar(ncid, 'domain_avg_kv_t', 'double', t_len);
netcdf.putAtt(ncid, domain_avg_kv_varid, 'description', ...
    sprintf('time series of domain-averaged Kv vertically averaged from %d m to bottom', vlevel_top));
netcdf.putAtt(ncid, domain_avg_kv_varid, 'units', 'm^2 s^-1');

box_avg_kv_varid = netcdf.defVar(ncid, 'box_avg_kv_t', 'double', t_len);
netcdf.putAtt(ncid, box_avg_kv_varid, 'description', ...
    sprintf('time series of box-averaged Kv vertically averaged from %d m to bottom', vlevel_top));
netcdf.putAtt(ncid, box_avg_kv_varid, 'units', 'm^2 s^-1');

netcdf.endDef(ncid);

% Write variables
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, domain_avg_kv_varid, domain_avg_kv_t);
netcdf.putVar(ncid, box_avg_kv_varid, box_avg_kv_t);

netcdf.close(ncid);

disp('Done saving vertically averaged Kv time series!')

