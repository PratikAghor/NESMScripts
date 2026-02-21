%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
% indxRange = 952:3877; % 
% indxRange = 3088:3095; % Jan 23, 2020
indxRange = 3024:3151; % Jan 15-30, 2020

nt0=indxRange(1);

rho_indxRange = 1691:3382;
rho_nt0 = rho_indxRange(1);
[~, rho_Nt] = size(rho_indxRange);

% Aghor's script
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
Kb = 1e-5; % background diffusivity
% Kmax = 5e-3;

%-------------------------------------------------------
sec = 'nozoom' % zoom or nozoom
%-------------------------------------------------------
%%
lonval = -63.15;
filename = sprintf('nesm_2019_2020_model_kv_t_vslice_const_lon_%.2f_nt_%d_%d.nc', lonval, indxRange(1), indxRange(Nt));

t_arr = ncread(filename, 'time');
lon_idx_arr = ncread(filename, 'lon_idx_arr');
lon_arr = ncread(filename, 'lon_arr');
kv_t_vslice = ncread(filename, 'kv_t_vslice_const_lon');
%%
%-------------------------------------------------------

%%
for nt = 1:Nt
    vort_t_nt = indxRange(nt) - nt0 + 1;

    for j = 1:length(lon_arr)
        % lon_idx = find(abs(lon_rho_vec-lon_arr(j))<1e-3); % find idx of lat in lon_rho_vec
        lon_idx = round(lon_idx_arr(j))
        figure1 = figure(j);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        Y = repmat(lat_rho_vec', NumLayers, 1);
        Z = (squeeze(zr(:, :, lon_idx)));
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        kv_vslice = squeeze(kv_t_vslice(nt, :, :));
        kv_vslice(kv_vslice<Kb) = Kb;
        % rho_vslice = squeeze(rho_pot_t_vslice(:, j, :));
    
        
        % rho_vslice = pagetranspose(rho_vslice);
    
        cmin = 1e-5; % min(min(vort_vslice)); 
        cmax = 1e-3; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        pcolor(ax1, Y, Z, kv_vslice); 
        shading interp;
        % freezeColors; hold on;
        set(ax1, 'ColorScale', 'log');
        colormap(ax1, cmocean('amp'));
        clim([cmin cmax]);
        
        % Colorbar with log ticks
        cb = colorbar;
        cb.Ticks = [1e-5, 1e-4, 1e-3];
        % cb.TickLabels = arrayfun(@(x) sprintf('%.e', x), cb.Ticks, 'UniformOutput', false);        
        cb.TickLabels = arrayfun(@(x) sprintf('$10^{%d}$', round(log10(x))), cb.Ticks, 'UniformOutput', false);
        set(cb, 'TickLabelInterpreter', 'latex');
        cb.FontSize = 14;
    
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, strcat('$$(k_{v}) \quad$$', date), 'interpreter','latex');
        hold on;
        
        % bathy
        plot(ax1, lat_rho_vec, squeeze(Z(1, :)), 'k-', 'LineWidth', 3);
        axis square
        % visual check for Symmetric instability (SI)
        if (strcmp(sec, 'zoom')==true)
            xlim(ax1,[38, 39]);
            ylim(ax1,[-5500, -2000]);
        elseif(strcmp(sec, 'nozoom')==true)
            % xlim([latlim]);
            ylim(ax1, [-5500, 0]);
        end
        %%%
        % 
        xticks([38 39 39.99]) 
        xticklabels({'38°N', '39°N','40°N'}) 
        yticks([-5000 -4000 -3000 -2000 -1000])     
        yticklabels({'5','4','3','2','1','0'});
        ylabel('Depth($\times 10^3$ m)', 'Interpreter', 'latex');
        %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        if (strcmp(sec, 'zoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lon_zoom/nesm_2019_2020_kv_vslice'];
        elseif(strcmp(sec, 'nozoom')==true)
            figname  = [plots_path, 'kv_plots/vslice/const_lon/nesm_2019_2020_kv_vslice'];
        end
        % 
        figname = strcat(figname, '_lon_', string(lon_arr(j)));
        figname = strcat(figname, '_nt_', string(indxRange(nt)));
        % 
        exportgraphics(figure1, strcat(figname, '.png')); % remove extra white space, 2022a and above

        close all;
    end
end
%------------------------------------

