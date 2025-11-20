clear all
% clc
%------------------------------------
% on prometheus
main_dir_path='/home/aghor/';
aghor_extras_path='/home/aghor/aghor/GT/GT_project/aghor_extras/'; % change according to the system
% on geometer
% main_dir_path='F:\prometheus_backup\Aghor\GT\GT_project\';
% aghor_extras_path = 'F:\prometheus_backup\Aghor\GT\GT_project\aghor_extras\';

addpath([aghor_extras_path, 'm_map/']);
plots_path = './';
%------------------------------------
% 
% lon = pagetranspose(ncread('NESM_grd.nc', 'lon_rho'));
% lat = pagetranspose(ncread('NESM_grd.nc', 'lat_rho'));
% depth = pagetranspose(ncread('NESM_grd.nc', 'h'));
% figname  = [plots_path, 'nesm_bathy'];
% 
% % 
% lon = pagetranspose(ncread('NESM_grd_5km.nc', 'lon_rho'));
% lat = pagetranspose(ncread('NESM_grd_5km.nc', 'lat_rho'));
% depth = pagetranspose(ncread('NESM_grd_5km.nc', 'h'));
% figname  = [plots_path, 'nesm_bathy'];

% lon = pagetranspose(ncread('NESM_grd_5km_to_1km.nc', 'lon_rho'));
% lat = pagetranspose(ncread('NESM_grd_5km_to_1km.nc', 'lat_rho'));
% depth = pagetranspose(ncread('NESM_grd_5km_to_1km.nc', 'h'));
% figname  = [plots_path, 'nesm_bathy_5km_to_1km'];

% lon = pagetranspose(ncread('SEAMOUNT_grd.nc', 'lon_rho'));
% lat = pagetranspose(ncread('SEAMOUNT_grd.nc', 'lat_rho'));
% depth = pagetranspose(ncread('SEAMOUNT_grd.nc', 'h'));
% figname  = [plots_path, 'seamount_bathy'];

lon = pagetranspose(ncread('SM3_grd.nc', 'lon_rho'));
lat = pagetranspose(ncread('SM3_grd.nc', 'lat_rho'));
depth = pagetranspose(ncread('SM3_grd.nc', 'h'));
figname  = [plots_path, 'sm3_bathy'];

% lon = pagetranspose(ncread('SM1_grd.nc', 'lon_rho'));
% lat = pagetranspose(ncread('SM1_grd.nc', 'lat_rho'));
% depth = pagetranspose(ncread('SM1_grd.nc', 'h'));
% figname  = [plots_path, 'sm1_bathy'];

lon_rho_vec = squeeze(lon(1, :));
lat_rho_vec = squeeze(lat(:, 1));
% 1km
lat_idx_min = 55;
lat_idx_max = 155;
lon_idx_min = 120;
lon_idx_max = 195;
% 5km
% lon_idx_min = 24;
% lon_idx_max = 40;
% lat_idx_min = 12;
% lat_idx_max = 33;

% Kelvin seamount idx
kelv_lat_idx_min = 120;
kelv_lat_idx_max = 180;
kelv_lon_idx_min = 60;
kelv_lon_idx_max = 120;

% Gosnold seamount idx
gosn_lat_idx_min = 10;
gosn_lat_idx_max = 130;
gosn_lon_idx_min = 195;
gosn_lon_idx_max = 280;

% %% Plotting
figure1 = figure();
[latlim,lonlim] = geoquadline(lat,lon);
m_proj('miller', 'long', lonlim,'lat', latlim);
zMin = 0;
zMax = 5000;
m_contourf(lon, lat, depth);
clim([zMin, zMax]);
m_grid('tickdir','in', ...
       'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
       'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
       'ytick',([38 39]), ... % latitude        
       'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;

hold on;
% box around Atlantis II


% box around Atlantis II
% plot_box(lat_rho_vec, lon_rho_vec, lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max, 'r');
% plot_box(lat_rho_vec, lon_rho_vec, kelv_lat_idx_min, kelv_lat_idx_max, kelv_lon_idx_min, kelv_lon_idx_max, 'b');
% plot_box(lat_rho_vec, lon_rho_vec, gosn_lat_idx_min, gosn_lat_idx_max, gosn_lon_idx_min, gosn_lon_idx_max, 'b');

% aghoredit

% title('NESM Bathymetry')
% title('Atlantis-II bathy')
h = colorbar;
% h.Label.String = "Depth (m)";
% h.Label.Rotation = 270;
% h.Label.VerticalAlignment = "bottom";
h.YTick = [1000 2000 3000 4000 5000 5500];
h.YTickLabel = {'1000', '2000', '3000', '4000', '5000', '5500'};
set(figure1, 'Visible', 'off'); % stop pop-ups

figname = strcat(figname, '.pdf');
exportgraphics(figure1, figname, 'ContentType', 'vector'); % remove extra white space, 2022a and above

%------------------------------------
%%
% % Create 3D plot
% figure2= figure('Visible', 'off'); % Set to 'on' to see the figure while running
% ax2 = axes(figure2);                  % Explicitly create axes
% s = surf(lon_rho_vec, lat_rho_vec, depth); % Flip depth to negative for elevation style
% % Appearance tweaks
% s.EdgeColor = 'none';                % Remove mesh lines
% s.FaceLighting = 'gouraud';          % Smooth lighting
% colormap(turbo);                     % Choose a vibrant colormap
% colorbar;
% clim([0 5000]);
% 
% % Add contour lines on top of the surface
% hold(ax2, 'on');
% [~, hContour] = contour3(ax2, lon, lat, depth, 30, 'k');  % 30 contour levels, black lines
% hContour.LineWidth = 0.5;
% 
% % Flip Z-axis so deeper depths are lower in 3D space
% set(ax2, 'ZDir', 'reverse');
% 
% xlabel('Longitude');
% ylabel('Latitude');
% zlabel('Depth (m)');
% title('3D Bathymetry: Atlantis II Seamount');
% 
% view(45, 30); % Adjust viewing angle
% axis tight;
% grid on;
% 
% % Optional: Export
% figname = strcat(plots_path, 'seamount_bathy_3d.pdf');
% exportgraphics(figure2, figname, 'ContentType', 'vector');
