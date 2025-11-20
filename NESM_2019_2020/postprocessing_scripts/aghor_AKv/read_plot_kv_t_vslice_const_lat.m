clear all; clc;
%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1
N = 1e-3; % Brunt-Vaisala frequency (order of magnitude, s^-1)
%%Loop to calculate averages if need (modify as needed)
indxRange = 952:3877; 
nt0=indxRange(1);
% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

rho_indxRange = 1691:3382;
rho_nt0 = rho_indxRange(1);
[~, rho_Nt] = size(rho_indxRange);

tskip = 1;
vlevel_skip = 50;
%-------------------------------------------------------
%-------------------------------------------------------
%%
kv_t_filename = strcat('nesm_2019_2020_model_kv_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

% t_arr = ncread(vort_t_filename, 't_arr');
lat_idx_arr = ncread(kv_t_filename, 'lat_idx_arr');
lat_arr = ncread(kv_t_filename, 'lat_arr');
kv_t_vslice = ncread(kv_t_filename, 'kv_t_vslice_const_lat');
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');
t_arr = ncread(kv_t_filename, 'time');
Kb = 1e-5; % background diffusivity
Kmax = 5e-3;
sec = 'nozoom' % zoom or nozoom
%%
%-------------------------------------------------------
% only to get the dates in the title!
%-------------------------------------------------------

% mkdir vort_plots;
%%
% 
%% plotting
for nt = 1:Nt
    vort_t_nt = indxRange(nt) - nt0 + 1;
    for i = 1:1 %length(lat_arr)
        % lat_idx = find(abs(lat_rho_vec-lat_arr(1))<1e-3); % find idx of lat in lat_rho_vec
    
        lat_idx = round(lat_idx_arr(i));
        figure1 = figure(i);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); % geoquadline requires Mapping Toolbox.
    
        X = (repmat(lon_rho_vec, NumLayers, 1));
        X = (X);
        Z = (squeeze(zr(:, lat_idx, :)));
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        kv_vslice = squeeze(kv_t_vslice(nt, :, i, :));
        kv_vslice(kv_vslice<Kb) = Kb;

        % kv_vslice(1, :) = 1e-4;
        % kv_vslice(kv_vslice>Kmax) = Kmax;
        cmin = 1e-5; % min(min(epv)); 
        cmax = 1e-3; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        % pcolor(ax1, Y(2:end-1, 2:end-1), Z(2:end-1, 2:end-1), epv_vslice(2:end-1, 2:end-1));
        pcolor(ax1, X, Z, kv_vslice); 
        shading interp;
        % freezeColors; hold on;
        set(ax1, 'ColorScale', 'log');
        colormap(ax1, cmocean('amp'));
        clim([cmin cmax]);
        
        % Colorbar with log ticks
        cb = colorbar;
        cb.Ticks = [1e-5, 1e-4, 1e-3];
        cb.TickLabels = arrayfun(@(x) sprintf('%.e', x), cb.Ticks, 'UniformOutput', false);        
        cb.FontSize = 14;

        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, date);
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        
        hold on;
        % Z = -flip(squeeze(zr(lat_idx,:,:)), 1);
        % LevelList1 = linspace(min(min(rho_vslice)), 36.49, 10);
        % LevelList2 = linspace(36.49, 36.8, 5);
        % LevelList3 = linspace(36.8, max(max(rho_vslice)), 10);
        % 
        % LevelList = [LevelList1, LevelList2, LevelList3];
        % contour(ax1, X, (Z), (rho_vslice), ...
        %     'EdgeColor', [0 0 0], 'LevelList', LevelList);
    
    
        % visual check for Symmetric instability (SI)
        if (strcmp(sec, 'zoom')==true)
            xlim([-63.5, -62.9]);
            ylim([-5000, -4000]);
        elseif(strcmp(sec, 'nozoom')==true)
            % xlim([-64., -61.]);
            ylim([-5000, 0]);
        end
        xticks([-64.99 -64 -63 -62 -61])  % longitude   
        xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
        yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
        axis square
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        if (strcmp(sec, 'zoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lat_zoom/nesm_2019_2020_model_kv_t_vslice'];
        elseif(strcmp(sec, 'nozoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lat/nesm_2019_2020_model_kv_t_vslice'];
        end
        % 
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(nt)));
        % vort_contour = strcat(vort_contour, '.pdf');
        exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        % set(figure1, 'PaperPositionMode', 'auto')
        % print(figure1,strcat(figname, '.png'),'-dpng','-r300');
        % exportgraphics(figure1, strcat(figname, '.eps'), 'Resolution',300); % remove extra white space, 2022a and above
        % export_fig(ax1, figname, '-eps','-transparent', '-r300'); % 
        close all;
    end
end
%------------------------------------

