%-------------------------------------------------------
%%
start_paths

%-------------------------------------------------------
%%
%%Loop to calculate averages if need (modify as needed)
% indxRange = 0:1690; % What time indices do you need?

indxRange = 952:3877; % entire year
% indxRange=3024:3031; % Jan 15, 2020
% indxRange=2032:2160; % Sep 13 - Sep 27, 2019
% indxRange=2232:2240; % Oct 08, 2019
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
% indxRange=3024:3151; % Jan 15 - Jan 30, 2020
% indxRange=3264:3392; % Feb 14 - Feb 29, 2020
% indxRange=3088:3097; % Jan 23, 2020
% indxRange=3136:3145; % Jan 29, 2020

uv_indxRange=indxRange;
vlevel = -4000;
% vlevel = -4000;
nt0=indxRange(1);
uv_nt0=uv_indxRange(1);

lon_min = -64.5; lon_max = -62.8;
lat_min = 38; lat_max = 39.5;
lonlim = [lon_min lon_max];
latlim = [lat_min lat_max];


% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
[~, uv_Nt] = size(uv_indxRange);

ke_arr = zeros(Nt, NumLayers); % get ke at each time at different depths
time_arr = zeros(Nt, 1);
%--------------------------------------------------------
tskip = 1;
vlevel_skip = 50;
filename = strcat('nesm_2019_2020_5km_kv_avg_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
% 
% t_arr = ncread(filename, 't_arr');
kv_avg =squeeze(ncread(filename, 'kv_avg_hslice'));

f0 = 0.909*1e-4; % Coriolis freq.
N0 = 1e-3; % Brunt Vaisala freq.
%%
%-------------------------------------------------------
%--------------------------------------------------------------------------
% vertical line profiles at const lat-lon vs z
% lat_idx_arr = ([90, 111, 135]);
lat_idx_arr = ([23]);

lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
%--------------------------------------------------------------------------
% lon_idx_arr = ([150, 158, 166, 172]);
lon_idx_arr = ([35]);
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
	lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j)); % const lon values to save vslice at
end
latsec = lat_rho_vec;
%--------------------------------------------------------

%-------------------------------------------------------
% mkdir vort_plots;
%% plotting
    %%% plotting
    figure1 = figure(1);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    % ax1 = axes('LineWidth', 3, 'FontSize', 24);

    m_proj('miller','long', lonlim,'lat', latlim);
    % h1 = m_pcolor(lon_rho, lat_rho, log10(R));
   

    hold on;
    kv = kv_avg;
    % u_hslice = squeeze(u_t(vort_t_nt, :, :));
    % v_hslice = squeeze(v_t(vort_t_nt, :, :));
    cmin = 1e-5;
    cmax = 1e-2; 
    h1 = m_image(lon_rho_vec, lat_rho_vec, kv);
    set(ax1, 'ColorScale', 'log');
    colormap(ax1, cmocean('amp')); clim([cmin cmax]);
    % cb = colorbar;
    % cb.Ticks = [1e-5, 1e-4, 1e-3, 1e-2];
    % cb.TickLabels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), cb.Ticks, 'UniformOutput', false);
    % set(cb, 'TickLabelInterpreter', 'latex');
    % cb.FontSize = 14;
    % 
    % xyskip = 15;
    % m_quiver(lon_rho(1:xyskip:end, 1:xyskip:end), lat_rho(1:xyskip:end, 1:xyskip:end), ...
    %  u_hslice(1:xyskip:end, 1:xyskip:end), v_hslice(1:xyskip:end, 1:xyskip:end), ...
    %  'color',[0 0 0]);
    % 
    % hold on;
    % % add a reference arrow using m_vec
    % [hpv5, htv5] = m_vec(100, -62, 40.25, 20, 0, 'k', 'key', '0.2 m/s');
    % % [hpv5, htv5] = m_vec(100, -61.5, 37.4, 20, 0, 'k', 'key', '0.2 m/s');
    % set(htv5,'FontSize',16);
    
    % [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [-3000 -3000], ...
    %     'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1 1 1], 'FaceAlpha', 0);
    [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
        'LineWidth', 1, 'EdgeColor', [0 0 0], 'FaceColor', [1 1 1], 'FaceAlpha', 0);
    % c2.FaceColor = [1 1 1]; c2.FaceAlpha = 0.3; % opacity
    % c2.EdgeColor = [0 0 0];
    % [h2, c2] = m_contour(lon_rho, lat_rho, -depth, [-3000 -3000], ...
    %        'LineWidth', 3, 'LineColor', [0 0 0]);
    
   %  m_grid('tickdir','in', ...
   % 'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
   % 'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
   % 'ytick',([38 39]), ... % latitude        
   % 'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
   
    m_grid('tickdir','in', ...
   'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
   'xticklabel',{}, ... % name longitude ticks as you want
   'ytick',([38 39]), ... % latitude        
   'yticklabel',{}); % name latitude ticks as you want;

    % --- Plot vertical section lines (constant lon) ---
    % for j = 1:length(lon_arr)
    %     m_line([lon_arr(j), lon_arr(j)], [min(lat_rho_vec(:)), max(lat_rho_vec(:))], ...
    %         'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
    % end
    
    % --- Plot horizontal section lines (constant lat) ---
    % for i = 1:length(lat_arr)
    %     m_line([min(lon_rho_vec(:)), max(lon_rho_vec(:))], [lat_arr(i), lat_arr(i)], ...
    %         'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
    % end

    % date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
    % date = date(1:11);
    % title(ax1, strcat('$$(k_v) \quad$$', date), 'interpreter','latex');
    title(ax1, strcat('Mean $$(k_v) \quad$$'), 'interpreter','latex');
    % ax1.TitleHorizontalAlignment = 'center'; % left makes it come to center
%%%
    set(figure1, 'Visible', 'off'); % stop pop-ups
    figname  = [plots_path, 'kv_plots/kv_avg_hslice/nesm_2019_2020_5km_kv'];
    
    figname = strcat(figname, '_vlevel_', string(vlevel));
    figname = strcat(figname, '_nt_', string(indxRange(1)), '_', string(indxRange(end)));
    figname = strcat(figname, '.png');
    exportgraphics(figure1, figname, 'ContentType', 'vector'); % remove extra white space, 2022a and above
    % 
    close all;
% 
%%

%------------------------------------

