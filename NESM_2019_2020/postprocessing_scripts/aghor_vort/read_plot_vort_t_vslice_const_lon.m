clear all; clc;
%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1
N = 1e-3; % Brunt-Vaisala frequency (order of magnitude, s^-1)
%%Loop to calculate averages if need (modify as needed)
% indxRange = 952:3877; % entire year

% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
indxRange=3024:3152; % Jan 15 - Jan 30, 2020
% indxRange=3264:3392; % Feb 14 - Feb 29, 2020

nt0=indxRange(1);
% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

% rho_indxRange = 1691:3382;
% rho_nt0 = rho_indxRange(1);
% [~, rho_Nt] = size(rho_indxRange);

tskip = 1;
% vlevel_skip = 50;
%-------------------------------------------------------
%-------------------------------------------------------
%%
vort_t_filename = strcat('nesm_2019_2020_vort_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

t_arr = ncread(vort_t_filename, 'time');
lon_idx_arr = ncread(vort_t_filename, 'lon_idx_arr');
lon_arr = ncread(vort_t_filename, 'lon_arr');
vort_t_vslice = ncread(vort_t_filename, 'vort_t_vslice_const_lon');

sec = 'nozoom' % zoom or nozoom
%%
%-------------------------------------------------------
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lon_nt_',string(indxRange(1)), '_', ...
%    string(indxRange(Nt)),'.nc');
% rho_pot_t_filename = strcat('nesm_2019_2020_rho_pot_t_vslice_const_lat_nt_',string(rho_indxRange(1)), '_', ...
%     string(rho_indxRange(rho_Nt)),'.nc');

% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_pot_t_vslice');
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');

%%
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------
% rho_pot_t_vslice = ncread(rho_pot_t_filename, 'rho_t_vslice');

% mkdir vort_plots;
%%
% 
%% plotting
for nt = 1:Nt
    vort_t_nt = indxRange(nt) - nt0 + 1;
    for j = 1:length(lon_arr)
        % lon_idx = find(abs(lon_rho_vec-lon_arr(j))<1e-3); % find idx of lat in lon_rho_vec
        lon_idx = round(lon_idx_arr(j))        
        
        figure1 = figure(j);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho); % geoquadline requires Mapping Toolbox.
    
        Y = repmat(lat_rho_vec', NumLayers, 1);
        Z = (squeeze(zr(:, :, lon_idx)));
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        vort_vslice = squeeze(vort_t_vslice(nt, :, :, j))./f0;
        vort_plus1 = vort_vslice + 1;
        % rho_pot_vslice = squeeze(rho_pot_t_vslice(nt, :, :, j));
        zMax = 2; % max(abs(vort_plus1(:)));
        zMin = -zMax;
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        % pcolor(ax1, Y(2:end-1, 2:end-1), Z(2:end-1, 2:end-1), vort_vslice(2:end-1, 2:end-1));
        pcolor(ax1, Y, Z, vort_vslice + 1); 
        shading interp;
        % freezeColors; hold on;
        set(ax1,'Color', [1 1 1])
        % set(ax1,'YDir','reverse')
        % colormap(ax1,b2r(zMin,zMax));  
        % colormap(ax1, whitejet); clim([zMin zMax]);
        % colormap(ax1, "jet"); 
        colormap(cmocean('curl')); cb = colorbar;
        cb.FontSize = 18;
        set(cb.Label, 'FontSize', 10)
        clim([zMin zMax]);
        colorbar;
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, strcat('$$(\zeta/f +1) \quad$$', date), 'interpreter','latex');
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        
        hold on;
        % LevelList1 = linspace(min(min(rho_pot_vslice)), 36.49, 20);
        % LevelList2 = linspace(36.49, 36.8, 10);
        % LevelList3 = linspace(36.8, max(max(rho_pot_vslice)), 20);

        % LevelList = [LevelList1, LevelList2, LevelList3];
        % contour(ax1, Y, (Z), (rho_pot_vslice), ...
        %    'EdgeColor', [0 0 0], 'LevelList', LevelList);
    
    
        % visual check for Symmetric instability (SI)
        if (strcmp(sec, 'zoom')==true)
            xlim([38, 39]);
            ylim([-5000, -4000]);
        elseif(strcmp(sec, 'nozoom')==true)
            % xlim([latlim]);
            ylim([-5000, 0]);
        end
        xticks([38 39 39.99])  % longitude   
        xticklabels({'38°N', '39°N','40°N'}) % name longitude ticks as you want
        yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
        axis square
        % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        if (strcmp(sec, 'zoom')==true)
            figname  = [plots_path, 'vort_plots/vslice/const_lon_zoom/nesm_2019_2020_vort_t_vslice'];
        elseif(strcmp(sec, 'nozoom')==true)
            figname  = [plots_path, 'vort_plots/vslice/const_lon/nesm_2019_2020_vort_t_vslice'];
        end
        % 
        figname = strcat(figname, '_lon_', string(lon_arr(j)));
        figname = strcat(figname, '_nt_', string(indxRange(nt)));
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

