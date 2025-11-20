function plot_vebf_hslice(case_name, vlevel, plots_path, lon_min, lon_max, lat_min, lat_max, ...
    lat_idx_arr, box, ...
    cbar, ref_arrow, cMax, xyskip, ...
    show_xlabel, show_ylabel, show_xticks, show_yticks)
%--------------------------------------------------------------------------
% Plot horizontal slice of VEBF with velocity arrows at a given depth.
% Inputs: same as plot_KmKe_hslice but for VEBF
%--------------------------------------------------------------------------
close all
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
aghor_extras_path = '../../../aghor_extras/';
addpath(fullfile(aghor_extras_path, 'export_fig'));
addpath(fullfile(aghor_extras_path, 'm_map'));
addpath(fullfile(aghor_extras_path, 'cmap_manual'));
%--------------------------------------------------------------------------
indxRange = 952:3877;
box = box;
[box_lat_idx_min, box_lat_idx_max, box_lon_idx_min, box_lon_idx_max] = box_idx_lims(box);

nt0 = indxRange(1); nt1 = indxRange(end);
vstr = string(vlevel);

VEBF_file = sprintf('%s_VEBF_2d_hslice_nt_%d_%d_vlevel_%s.nc', case_name, nt0, nt1, vstr);
uv_file   = sprintf('../uv_annual_means/%s_annual_mean_uv_hslice_vlevel_%s.nc', case_name, vstr);


VEBF2D   = ncread(fullfile(plots_path, VEBF_file), 'VEBF2D');
u_hslice = ncread(uv_file, 'u_hslice');
v_hslice = ncread(uv_file, 'v_hslice');
lon_rho  = ncread(uv_file, 'lon_rho');
lat_rho  = ncread(uv_file, 'lat_rho');

lon_rho_vec = squeeze(lon_rho(1,:));
lat_rho_vec = squeeze(lat_rho(:,1));

lonlim = [lon_rho_vec(1), lon_rho_vec(end)];
latlim = [lat_rho_vec(1), lat_rho_vec(end)];

if nargin >= 4 && ~isempty(lon_min) && ~isempty(lon_max)
    lonlim = [lon_min lon_max];
end
if nargin >= 6 && ~isempty(lat_min) && ~isempty(lat_max)
    latlim = [lat_min lat_max];
end

figure1 = figure('Visible', 'off');
axes('Parent', figure1, 'YMinorTick', 'on', 'LineWidth', 2, 'FontSize', 20);

% Use m_proj with view limits only
m_proj('miller', 'long', lonlim, 'lat', latlim);
hold on;

cMin = -cMax;

[lon_mesh, lat_mesh] = meshgrid(lon_rho_vec, lat_rho_vec);
m_pcolor(lon_mesh, lat_mesh, VEBF2D);
shading interp;
colormap(cmocean('curl'));
clim([cMin cMax]);
axis tight

if cbar
    cb = colorbar;
    cb.FontSize = 16;
else
    colorbar('Visible', 'off');
end

xyskip = max(1,xyskip);

m_quiver(lon_rho(1:xyskip:end, 1:xyskip:end), ...
         lat_rho(1:xyskip:end, 1:xyskip:end), ...
         u_hslice(1:xyskip:end, 1:xyskip:end), ...
         v_hslice(1:xyskip:end, 1:xyskip:end), 'k');

%--------------------------------------------------------------------------
% horizontal lines at given lat_idx_arr

lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end

% Plot horizontal section lines (constant lat)
% for i = 1:length(lat_arr)
%     m_line([min(lon_rho_vec(:)), max(lon_rho_vec(:))], [lat_arr(i), lat_arr(i)], ...
%        'linestyle', '--', 'Color', 'blue', 'linewidth', 3);
% end
%--------------------------------------------------------------------------
plot_box(lat_rho_vec, lon_rho_vec, box_lat_idx_min, box_lat_idx_max, box_lon_idx_min, box_lon_idx_max, '-b');
%--------------------------------------------------------------------------
if ref_arrow
    lon_pos = lonlim(2) - 0.4;
    lat_pos = latlim(2) + 0.12;
    [~, htv] = m_vec(100, lon_pos, lat_pos, 20, 0, 'k', 'key', '0.2 m/s');
    set(htv, 'FontSize', 16);
end

if show_xticks
    % Generate ticks that are actually within plot bounds
    lon_range = lonlim(2) - lonlim(1);
    lon_step = 1.0;  % 1 degree steps for wide ranges
    
    % Find ticks within the actual bounds
    first_tick = ceil(lonlim(1) / lon_step) * lon_step;
    last_tick = floor(lonlim(2) / lon_step) * lon_step;
    xticks_vals = first_tick:lon_step:last_tick;
    
    % Generate labels
    xticklabels_vals = arrayfun(@(x) sprintf('%.0f°W', abs(x)), xticks_vals, 'UniformOutput', false);
else
    xticks_vals = [];
    xticklabels_vals = {};
end

if show_yticks
    % Generate ticks that are actually within plot bounds
    lat_range = latlim(2) - latlim(1);
    lat_step = 1.0;  % 1 degree steps for wide ranges
    
    % Find ticks within the actual bounds
    first_tick = ceil(latlim(1) / lat_step) * lat_step;
    last_tick = floor(latlim(2) / lat_step) * lat_step;
    yticks_vals = first_tick:lat_step:last_tick;
    
    % Generate labels
    yticklabels_vals = arrayfun(@(y) sprintf('%.0f°N', y), yticks_vals, 'UniformOutput', false);
else
    yticks_vals = [];
    yticklabels_vals = {};
end

m_grid('tickdir','in', ...
    'xtick', xticks_vals, 'xticklabel', xticklabels_vals, ...
    'ytick', yticks_vals, 'yticklabel', yticklabels_vals, ...
    'xminortick','off', 'yminortick','off');

if show_xlabel
    xlabel('');
else
    xlabel('');
end

if show_ylabel
    ylabel('');
else
    ylabel('');
end

% Set figure size to avoid PDF page size error
set(figure1, 'Units', 'inches', 'Position', [0, 0, 5, 5]);
figname = sprintf('%s_VEBF_2d_hslice_nt_%d_%d_vlevel_%s.pdf', ...
    case_name, nt0, nt1, vstr);

exportgraphics(figure1, fullfile(plots_path, figname), ...
    'ContentType', 'vector', 'Resolution', 300);

close(figure1);
fprintf('Saved: %s\n', fullfile(plots_path, figname));
%--------------------------------------------------------------------------
end
