% 
% Compute the kinetic energy transfer: horizontal Reynolds stress (HRS)
% - Pratik Aghor, modified from get_KmKe.m
% HRS= -[ <up up> dubar/dx + <up vp>  dubar/dy  ....
%          <vp up> dvbar/dx + <vp vp>  dvbar/dy ]
%     =  -<up(up.grad(ubar))+vp(up.grad(vbar))>
%
% Advection operators are used. Much easier in sigma coordinates.
%
% This is the method used in
% Djakouré, S., P. Penven, B. Bourlès, J. Veitch and V. Koné, 
% Coastally trapped eddies in the north of the Gulf of Guinea, 2014, J. Geophys. Res,
% DOI: 10.1002/2014JC010243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Assume mean and hrs_3d with respect to that mean has been calculated and saved
% only take a horizontal section of the saved hrs_3d at a given depth, 
% and save HRS2D 
 
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % What time indices do you need?
% indxRange = 2288:2416; % Oct. 15-30, 2019
% indxRange = 3024:3152; % Jan. 15-30, 2020

vlevel = -3000;
nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);

KmKe_2d_indxRange = indxRange;
[~, mean_Nt] = size(KmKe_2d_indxRange);

uvw_indxRange = 952:3877;
    [~, uvw_Nt] = size(uvw_indxRange);

VEBF_2d_file = strcat('sm3_2019_2020_VEBF_2d_hslice_nt_', string(KmKe_2d_indxRange(1)), '_', string(KmKe_2d_indxRange(mean_Nt)), ...
    '_vlevel_', string(vlevel), '.nc');

uvw_file = strcat('../aghor_kmke/sm3_2019_2020_uvw_annual_mean_3d_nt_', string(uvw_indxRange(1)), '_', string(uvw_indxRange(uvw_Nt)), '.nc');
%-------------------------------------------------------
VEBF2D = ncread(VEBF_2d_file, 'VEBF2D');
ub = u2rho_3d(ncread(uvw_file, 'ub'));
vb = v2rho_3d(ncread(uvw_file, 'vb'));
u_hslice = vinterp(ub, zr, vlevel);
v_hslice = vinterp(vb, zr, vlevel);
%-------------------------------------------------------

%--------------------------------------------------------------------------
%-------------------------------------------------------

%--------------------------------------------------------------------------

%% plot the time average, real KmKe2d, now already done in the modified save_KmKe_2d.m script
% KmKe2d_avg = zeros(Ny, Nx);
% for nt=1:Nt
%     KmKe2d_avg = KmKe2d_avg + squeeze(KmKe_t(nt, :, :));
% end
% KmKe2d_avg = KmKe2d_avg./Nt;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% plot VEBF hslice
% box_arr = ["N", "NE", "W", "SW", "S", "SE", "E", "NW"]
% box_arr = ["big_N", "big_S"]
box_arr = []

figure1 = figure(3);
% [latlim, lonlim] = geoquadline(lat_rho, lon_rho); geoquadline requires
% Mapping toolbox
% Create axes
ax1 = axes('Parent', figure1, 'YMinorTick','on',...
    'LineWidth',3,...
    'FontSize',24);

m_proj('miller','long', lonlim,'lat', latlim);


hold on;
zMin = -1e-7; % min(min(vort)); 
zMax = -zMin;
h1 = m_image(lon_rho_vec, lat_rho_vec, VEBF2D);
colormap('parula'); 
cb = colorbar; 
cb.FontSize = 20;
clim([zMin zMax]);
% colormap(ax1,b2r(zMin,zMax));  colorbar;
% colormap(ax1, "jet"); clim([zMin zMax]); colorbar;
% colormap(ax1, whitejet); clim([zMin zMax]); colorbar;
% freezeColors; hold on;

[h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
    'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1  1 1], 'FaceAlpha', 1);
% c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
% c2.EdgeColor = [0 0 0];

xyskip = 15;
m_quiver(lon_rho(1:xyskip:end, 1:xyskip:end), lat_rho(1:xyskip:end, 1:xyskip:end), ...
 u_hslice(1:xyskip:end, 1:xyskip:end), v_hslice(1:xyskip:end, 1:xyskip:end), ...
 'color',[0 0 0]);

hold on;
% add a reference arrow using m_vec
[hpv5, htv5] = m_vec(100, -62, 40.25, 20, 0, 'k', 'key', '0.2 m/s');
% [hpv5, htv5] = m_vec(100, -61.5, 37.4, 20, 0, 'k', 'key', '0.2 m/s');
set(htv5,'FontSize',16);

m_grid('tickdir','in', ...
'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
'ytick',([38 39]), ... % latitude        
'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;

% title(ax1, string(datetime(2017,12,31,0,0,t_arr(vort_t_nt))));
% ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
hold on;
% title(ax1, string(date));
% plot boxes over which horizontal averages will be taken
% for i = 1:length(box_arr)
%     box = string(box_arr(i))
%     [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box);
% 
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_min)], ...
%     'r', 'LineWidth', 3);
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_max) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% 
% m_plot([lon_rho_vec(lon_idx_min) lon_rho_vec(lon_idx_min)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% m_plot([lon_rho_vec(lon_idx_max) lon_rho_vec(lon_idx_max)],[lat_rho_vec(lat_idx_min) lat_rho_vec(lat_idx_max)], ...
%     'r', 'LineWidth', 3);
% 
% hold on;
% end
%%%
set(figure1, 'Visible', 'off'); % stop pop-ups
figname  = [plots_path, 'rs_plots/vebf/sm3_2019_2020_vebf_2d_hslice'];

figname = strcat(figname, '_vlevel_', string(vlevel));
figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(Nt)));
exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
% exportgraphics(figure1,strcat(figname, '.eps'))
close all;
%--------------------------------------------------------------------------