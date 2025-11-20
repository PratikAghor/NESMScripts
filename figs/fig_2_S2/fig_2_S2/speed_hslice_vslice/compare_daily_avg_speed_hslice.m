%-------------------------------------------------------
%%
start_paths;
if ~exist('daily_avg_speed_hslice','dir'); mkdir('daily_avg_speed_hslice'); end
outdir = 'daily_avg_speed_hslice';
%-------------------------------------------------------
indxRange = 952:3877; % entire year
% indxRange = 2912:3023; % January 1-14, 2020

uv_indxRange=indxRange;
% vlevel = -500;
vlevel = -4000;

cMin = 0;
if(vlevel == -500 & indxRange(1) == 952)
    cMax = 0.2;
elseif(vlevel == -500 & indxRange(1) ~= 952)
    cMax = 0.3;
elseif(vlevel == -4000)
    cMax = 0.05;
end

nt0=indxRange(1);
uv_nt0=uv_indxRange(1);

% Aghor's script
[~, Nt] = size(indxRange);
[~, uv_Nt] = size(uv_indxRange);

%--------------------------------------------------------------------------

%--------------------------------------------------------
% Model prefixes and their dirs
models = { ...
  struct('prefix','nesm_1km', 'filename_prefix', 'nesm_2019_2020', ...
  'speed_path','.'), ...
  % struct('prefix','nesm_5km', 'filename_prefix', 'nesm_2019_2020_5km', ...
  % 'speed_path','.'), ...
  % struct('prefix','sm3_1km',  'filename_prefix', 'sm3_2019_2020', ...
  % 'speed_path','.'), ...
};

% containers
avg_speed_hslice_all = cell(numel(models),1);
avg_ke_hslice_all = cell(numel(models),1);

avg_u_hslice_all = cell(numel(models),1);
avg_v_hslice_all = cell(numel(models),1);

labels      = cell(numel(models),1);

for m = 1:numel(models)
    prefix   = models{m}.prefix;
    filename_prefix   = models{m}.filename_prefix;
    speed_path   = models{m}.speed_path;

    c_t_filename = sprintf('%s/%s_avg_speed_t_hslice_vlevel_%d_nt_%d_%d.nc', ...
        speed_path, filename_prefix, vlevel, indxRange(1), indxRange(Nt));
    
    avg_speed_hslice =squeeze(ncread(c_t_filename, 'avg_speed_t'));
    avg_ke_hslice =squeeze(ncread(c_t_filename, 'avg_ke_t'));
    avg_u_hslice =squeeze(ncread(c_t_filename, 'avg_u_t'));
    avg_v_hslice =squeeze(ncread(c_t_filename, 'avg_v_t'));
    
    avg_speed_hslice_all{m} = avg_speed_hslice;
    avg_ke_hslice_all{m} = avg_ke_hslice;
    avg_u_hslice_all{m} = avg_u_hslice;
    avg_v_hslice_all{m} = avg_v_hslice;

    labels{m}      = upper(strrep(prefix,'_',' '));
end
%%
%-------------------------------------------------------
depth = nesm_1km_depth;
%--------------------------------------------------------------------------
%% KE
for m = 1:numel(models)
    prefix = models{m}.prefix;
    

    if strcmp(prefix, 'nesm_1km')
        lon_rho = lon_rho_1km;
        lat_rho = lat_rho_1km;
        depth = nesm_1km_depth;
        xyskip = 15;
        lat_idx_arr = ([111]);
        lon_idx_arr = ([159])


    elseif strcmp(prefix, 'nesm_5km')
        lon_rho = lon_rho_5km;
        lat_rho = lat_rho_5km;
        depth = nesm_5km_depth;
        xyskip = 3;
        lat_idx_arr = ([23]);
        lon_idx_arr = ([33])

    elseif strcmp(prefix, 'sm3_1km')
        lon_rho = sm3_1km_lon_rho;
        lat_rho = sm3_1km_lat_rho;
        depth = sm3_1km_depth;
        xyskip = 15;
        lat_idx_arr = ([111]);
        lon_idx_arr = ([159]);
    end
    lon_rho_vec= squeeze(lon_rho(1, :));
    lat_rho_vec = squeeze(lat_rho(:, 1));
    %--------------------------------------------------------------------------
    % vertical line at const lat-lon vs z

    lat_arr = zeros(length(lat_idx_arr), 1); % const lat values to save vslice at
    for i = 1:length(lat_idx_arr)
            lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
    end

    % horizontal line at const lat-lon vs z

    lon_arr = zeros(length(lon_idx_arr), 1); % const lat values to save vslice at
    for j = 1:length(lon_idx_arr)
            lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));
    end

    %--------------------------------------------------------------------------

    lonlim = [lon_rho_vec(1), lon_rho_vec(end)];
    % lonlim = ([-64 -62]);
    latlim = [lat_rho_vec(1), lat_rho_vec(end)];
    
    % plot avg speed hslice
    avg_ke_hslice = avg_ke_hslice_all{m};
    avg_u_hslice = avg_u_hslice_all{m};
    avg_v_hslice = avg_v_hslice_all{m};
    
    % Plot ke hslice
    fig = figure('Visible','off');
    set(fig, 'Units', 'inches', 'Position', [0, 0, 5, 5]);
    ax1 = axes('Parent', fig, 'YMinorTick','on', 'LineWidth',3, 'FontSize',24);

    m_proj('miller','long', lonlim,'lat', latlim); hold(ax1,'on');

    % contourf(lon_rho, lat_rho, avg_speed_hslice, 50, 'LineColor','none');
    m_image(lon_rho_vec, lat_rho_vec, avg_ke_hslice);
    % colormap(ax1, cmocean('speed'));
    colormap(ax1, parula(8));
    cMin = 0;
    % if(vlevel == -500)
    %     cMax = 0.2;
    % elseif(vlevel == -4000)
    %     cMax = 0.1
    % end

    clim(ax1, [cMin cMax]);
    arrow_scale_factor = 1;
    if(strcmp(prefix, 'nesm_1km'))
        % Colorbar
        % colorbar('peer', ax1);
        %-----------------
        cb = colorbar(ax1);
        cb.FontSize = 12;
        cb.Ticks = [cMin, cMax];
        if(vlevel == -500)
            cb.Ruler.Exponent = -1; % 10^{-1}
            cb.Ruler.TickLabelFormat = '%.0f';
        elseif(vlevel == -4000)
            cb.Ruler.Exponent = -2; % 10^{-2}
            cb.Ruler.TickLabelFormat = '%.0f';
        end 
    end
    hold on;
    if(strcmp(prefix, 'nesm_1km'))
        % xlabel('Longitude'); ylabel('Latitude');
        xtick_vals = ([-64.99 -64 -63 -62]);
        xticklabels = {'65°W','64°W','63°W','62°W'};
        ytick_vals = ([38 39]);
        yticklabels = {'38°N','39°N'};
    else
        xtick_vals = ([]);
        xticklabels = {};
        ytick_vals = ([]);
        yticklabels = {};
    end
    set(ax1, 'FontSize', 18);           
    % Add velocity quivers
    m_quiver(lon_rho(1:xyskip:end,1:xyskip:end), ...
           lat_rho(1:xyskip:end,1:xyskip:end), ...
           avg_u_hslice(1:xyskip:end,1:xyskip:end), ...
           avg_v_hslice(1:xyskip:end,1:xyskip:end), ...
           arrow_scale_factor, 'color','k');
    
    
    hold on
    m_grid('tickdir','in', ...
            'xtick', xtick_vals, ...
            'xticklabel', xticklabels, ...
            'ytick', ytick_vals, ...
            'yticklabel', yticklabels, ...
            'aspect', 1);

    % Bathy contour
    m_contourf(lon_rho, lat_rho, -depth, [-4000 -4000], ...
        'LineWidth',3, 'EdgeColor', [0 0 0], 'FaceColor', [1 1 1], 'FaceAlpha',0);
    
        %  Plot section lines (constant lat, lon)
        for i = 1:length(lat_arr)
            m_line([min(lon_rho_vec(:)), max(lon_rho_vec(:))], [lat_arr(i), lat_arr(i)], ...
                'color', 'w', 'linestyle', '--', 'linewidth', 1.5);
        end
        
        %  Plot vertical section lines (constant lon) 
        for j = 1:length(lon_arr)
            m_line([lon_arr(j), lon_arr(j)], [min(lat_rho_vec(:)), max(lat_rho_vec(:))], ...
                'color', 'w', 'linestyle', '--', 'linewidth', 1.5);
        end


    %-----------------
    % ref arrow 

    if(strcmp(prefix, 'nesm_1km'))
        % Position the reference arrow
        % lon_pos = lonlim(2) - 0.8;
        % lat_pos = latlim(2) + 0.2;
        lon_pos = lonlim(2) - 0.8;
        lat_pos = latlim(2) - 0.1;
        u_ref = 0.1; % m/s
        v_ref = 0;
       
       patch_width = 0.8;
       patch_height = 0.08;
       m_patch([lon_pos-0.05 lon_pos+patch_width lon_pos+patch_width lon_pos-0.05], ...
            [lat_pos-0.18 lat_pos-0.18 lat_pos+patch_height lat_pos+patch_height], ...
            'w', 'FaceAlpha', 0.75, 'EdgeColor', 'k', 'LineWidth', 1);

       hq = m_quiver(lon_pos, lat_pos, u_ref, v_ref, arrow_scale_factor, ...
        'color', 'k', 'LineWidth', 1, 'AutoScale', 'off', 'Clipping', 'off', ...
        'MaxHeadSize', 3);
        uistack(hq, 'top');
    
       m_text(lon_pos, lat_pos - 0.08, '0.1 m/s', ...
            'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    end
    %-----------------   
    figname = sprintf('%s/%s_avg_ke_t_hslice_vlevel_%d_nt_%d_%d.png', ...
        outdir, prefix, vlevel, indxRange(1), indxRange(Nt));

    exportgraphics(fig, figname, 'ContentType','vector');
    close(fig);
end