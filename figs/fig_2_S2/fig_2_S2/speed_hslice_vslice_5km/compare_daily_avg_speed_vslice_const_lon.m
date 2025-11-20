%-------------------------------------------------------
%%
start_paths;
if ~exist('daily_avg_speed_hslice','dir'); mkdir('daily_avg_speed_hslice'); end
outdir = 'daily_avg_speed_hslice';
%-------------------------------------------------------
indxRange = 952:3877; % entire year
cMin = 0; cMax = 0.2;
% indxRange = 2912:3023; % January 1-14, 2020
% cMin = 0; cMax = 0.3;

nt0=indxRange(1);


% Aghor's script
[~, Nt] = size(indxRange);
hline_depth = -500;
hline_depth_2 = -4000;
%--------------------------------------------------------
% Model prefixes and their dirs
models = { ...
  % struct('prefix','nesm_1km', 'filename_prefix', 'nesm_2019_2020', ...
  % 'speed_path','.', 'lonval', -63.15), ...
  struct('prefix','nesm_5km', 'filename_prefix', 'nesm_2019_2020_5km', ...
  'speed_path','.', 'latval', 38.50, 'lonval', -63.14), ...
  % struct('prefix','sm3_1km',  'filename_prefix', 'sm3_2019_2020', ...
  % 'speed_path','.', 'latval', 38.50) ...
};

% containers
avg_speed_vslice_all = cell(numel(models),1);
avg_ke_vslice_all = cell(numel(models),1);
avg_Z_vslice_all = cell(numel(models),1);

labels      = cell(numel(models),1);

for m = 1:numel(models)
    prefix   = models{m}.prefix;
    filename_prefix   = models{m}.filename_prefix;
    speed_path   = models{m}.speed_path;
    lonval = models{m}.lonval;

    c_t_filename = sprintf('%s/%s_avg_speed_t_vslice_const_lon_%.2f_nt_%d_%d.nc', ...
        speed_path, filename_prefix, lonval, indxRange(1), indxRange(Nt));
    disp(c_t_filename)
    avg_speed_vslice =squeeze(ncread(c_t_filename, 'avg_speed_t_const_lon'));
    avg_ke_vslice =squeeze(ncread(c_t_filename, 'avg_ke_t_const_lon'));
    % disp(size(avg_speed_vslice))

    Z =squeeze(ncread(c_t_filename, 'Z_slice'));
    
    avg_speed_vslice_all{m} = avg_speed_vslice;
    avg_ke_vslice_all{m} = avg_ke_vslice;
    avg_Z_vslice_all{m} = Z;

    labels{m}      = upper(strrep(prefix,'_',' '));
end
%%
%-------------------------------------------------------
depth = nesm_1km_depth;
%--------------------------------------------------------------------------
%% KE
for m = 1:numel(models)
    prefix = models{m}.prefix;
    avg_ke_vslice = avg_ke_vslice_all{m};
    Z = avg_Z_vslice_all{m};
    disp(size(avg_speed_vslice))
    
    if strcmp(prefix, 'nesm_1km')
        lon_rho = lon_rho_1km;
        lat_rho = lat_rho_1km;
        depth = nesm_1km_depth;
        xyskip = 15;

    elseif strcmp(prefix, 'nesm_5km')
        lon_rho = lon_rho_5km;
        lat_rho = lat_rho_5km;
        depth = nesm_5km_depth;
        xyskip = 3;

    elseif strcmp(prefix, 'sm3_1km')
        lon_rho = sm3_1km_lon_rho;
        lat_rho = sm3_1km_lat_rho;
        depth = sm3_1km_depth;
        xyskip = 15;
    end
    lon_rho_vec= squeeze(lon_rho(1, :));
    % disp(length(lon_rho_vec))
    lat_rho_vec = squeeze(lat_rho(:, 1));
    

    fig =  figure('Visible','off');
    ax = axes('YMinorTick','on','LineWidth',2,'FontSize',18); hold on;

    Y = repmat(lat_rho_vec', size(Z,1), 1);
    pcolor(ax, Y, Z, avg_ke_vslice); 
    shading interp;
    % colormap(cmocean('speed'));
    colormap(ax, parula(8));
    
    clim([cMin cMax]);
    % xlim([min(lon_rho_vec) max(lon_rho_vec)]);
    % xlim([-64 -62]);
    % pbaspect([1 1.5 1]);
    % pbaspect([1 1 1]);

    % ticks and labels
    if strcmp(prefix, 'nesm_5km')
        xticks([38, 39]);
        xticklabels({'38°N', '39°N'});
        yticks([-5000 -4000 -3000 -2000 -1000 0]);
        yticklabels({'5','4','3','2','1','0'});
        ylabel('Depth($\times 10^3$ m)', 'Interpreter', 'latex'); 
    else
        xticks([-65 -64 -63 -62 -61]); xticklabels({});
        yticks([-5000 -4000 -3000 -2000 -1000 0]); yticklabels({});
    end
    set(ax, 'FontSize', 18); 
    % colorbar
    if(strcmp(prefix, 'nesm_5km'))
        cb = colorbar(ax);
        cb.FontSize = 12;
        cb.Ticks = [cMin, cMax];
        cb.Ruler.Exponent = -1; % 10^{-1}
        cb.Ruler.TickLabelFormat = '%.0f';
    end
    
    % layout and save
    set(ax,'Color',[1 1 1]);
    
    % hline
    plot([lat_rho_vec(1) lat_rho_vec(end)], [hline_depth hline_depth], 'w--', 'LineWidth', 1.5);
    plot([lat_rho_vec(1) lat_rho_vec(end)], [hline_depth_2 hline_depth_2], 'w--', 'LineWidth', 1.5);

    % bathy
    plot(ax, lat_rho_vec, squeeze(Z(1, :)), 'k-', 'LineWidth', 3);
    axis tight;
    ylim(ax, [-5500 0]);

    figname = sprintf('%s/%s_avg_ke_vslice_lon_%.2f_nt_%d_%d.png', ...
        outdir, prefix, models{m}.lonval, indxRange(1), indxRange(end));
    exportgraphics(fig, figname, 'Resolution', 300);
    close(fig);
end
