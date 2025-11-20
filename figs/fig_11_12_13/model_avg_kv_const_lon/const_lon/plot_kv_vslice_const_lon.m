function plot_kv_vslice_const_lon(case_name, lonval, indxRange, plots_path, ...
    lat_min, lat_max, z_min, z_max, cbar, lat_idx, ...
    show_xlabels, show_ylabels)

%--------------------------------------------------------------------------
% plot_kv_vslice_const_lon
%
% Plot vertical slice of time-averaged Kv at constant longitude.
%
% Inputs:
%   case_name     : e.g., 'nesm_2019_2020'
%   lonval        : -63.25
%   indxRange     : time indices used in averaging, for filename
%   plots_path    : folder to save plots
%   lat_min/max   : x-axis limits (in degrees North)
%   z_min/max     : y-axis limits (in meters, negative downward)
%   cbar          : true/false, show colorbar
%   lat_idx       : optional index to draw vertical line at a latitude (set [] to skip)
%   show_xlabels  : true/false
%   show_ylabels  : true/false
%--------------------------------------------------------------------------

%--------------------------------------------------
% Color scale and background threshold
Kb   = 1e-5;
cmin = 1e-5;
cmax = 1e-3;

%--------------------------------------------------
% Construct filename and read data
file = sprintf('%s_avg_kv_t_vslice_const_lon_%.2f_nt_%d_%d.nc', ...
               case_name, lonval, indxRange(1), indxRange(end));

kv     = ncread(file, 'avg_kv_t_const_lon');  % [Nz x Ny]
Z      = ncread(file, 'Z_slice');             % [Nz x Ny]
lat_rho_vec    = ncread(file, 'lat_rho_vec');         % [Ny x 1]
lon    = ncread(file, 'lon_arr');             % [1 x 1]
[Nz, Ny] = size(Z);

kv(kv < Kb) = Kb;
% kv = kv;         % Now [Nz x Ny] to match Z
% Z = Z;           % [Nz x Ny] to match kv
Y = repmat(lat_rho_vec', Nz, 1);
size(Y)
size(Z)
size(kv)
%--------------------------------------------------------------------------
% aghor_extras_path = '../../../aghor_extras/';
% addpath(fullfile(aghor_extras_path, 'export_fig'));
% addpath(fullfile(aghor_extras_path, 'cmap_manual/'));
%--------------------------------------------------------------------------
%--------------------------------------------------
% Plot
figure('Visible', 'off');
ax = axes('FontSize', 18, 'LineWidth', 1.5);

pcolor(Y, Z, kv); 
shading interp;
set(ax, 'YDir', 'normal', 'ColorScale', 'log');
colormap(ax, cmocean('amp'));
clim([cmin cmax]);

%--------------------------------------------------
% Optional vertical line at lat_idx
if ~isempty(lat_idx)
    hold on;
    plot([lat_rho_vec(lat_idx), lat_rho_vec(lat_idx)], vline_depth, 'k--', 'LineWidth', 2);
end

%--------------------------------------------------
% Axis limits and ticks
xlim([lat_min, lat_max]);
ylim([z_min, z_max]);



if show_xlabels
    xticks([38 39 40])
    xticklabels({'38°N', '39°N', '40°N'})
    xlabel('');
else
    xticks([38 39 40]);
    ax.XTickLabel = [];
end

if show_ylabels
    ytick_vals = [-5000 -4000 -3000 -2000 -1000 0];
    yticks(ytick_vals);
    yticklabels({'5', '4', '3', '2', '1', '0'})
    ylabel('Depth($\times 10^3$ m)', 'Interpreter', 'latex', 'FontSize', 24);        % ylabel('Depth (m)');
else    
    ax.YTickLabel = [];
end
ax.FontSize = 18;
hold on;
% bathy
plot(ax, lat_rho_vec, squeeze(Z(1, :)), 'k-', 'LineWidth', 3); hold on;
% axis square
% Aspect ratio width:height = 1:1.5 (tall)
pbaspect([1 1.5 1]);

%--------------------------------------------------
% Colorbar
if cbar
    cb = colorbar(ax);
    cb.Ticks = [1e-5 1e-4 1e-3];
    cb.TickLabels = {'10^{-5}', '10^{-4}', '10^{-3}'};
    cb.FontSize = 14;
    % ylabel(cb, '$$k_v\ [\mathrm{m}^2/\mathrm{s}]$$', 'Interpreter', 'latex');
end

%--------------------------------------------------
% Title
% title(ax, sprintf('%s | lon = %.2f°W', strrep(case_name, '_', '\_'), lon), ...
%       'Interpreter', 'none');

%--------------------------------------------------
% Save
if ~exist(plots_path, 'dir'), mkdir(plots_path); end
figname = sprintf('%s%s_kv_vslice_nt_%d_%d.pdf', plots_path, case_name, ...
    indxRange(1), indxRange(end));
exportgraphics(gcf, figname, 'ContentType', 'vector');
close(gcf);

end

