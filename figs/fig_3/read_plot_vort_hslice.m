%-------------------------------------------------------
%%
start_paths;
%-------------------------------------------------------
%%
indxRange=3112:3119; % Jan 26, 2020
% config = 'nesm_2019_2020';
config = 'nesm_2019_2020_5km';
% config = 'sm3_2019_2020';

uv_indxRange=indxRange;
vlevel = -4000;
% vlevel = -4000;
nt0=indxRange(1);
uv_nt0=uv_indxRange(1);

% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
[~, uv_Nt] = size(uv_indxRange);

time_arr = zeros(Nt, 1);
%--------------------------------------------------------
tskip = 1;
vlevel_skip = 50;
vort_t_filename = strcat(config,'_vort_t_hslice_vlevel_', string(vlevel), '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
u_t_filename = strcat(config,'_u_t_hslice_vlevel_', string(vlevel), '_nt_',string(uv_indxRange(1)), '_', ...
    string(uv_indxRange(uv_Nt)),'.nc');

v_t_filename = strcat(config,'_v_t_hslice_vlevel_', string(vlevel), '_nt_',string(uv_indxRange(1)), '_', ...
    string(uv_indxRange(uv_Nt)),'.nc');
% 
t_arr = ncread(vort_t_filename, 't_arr');
vort_t =squeeze(ncread(vort_t_filename, 'vort_t_hslice'));
u_t = squeeze(ncread(u_t_filename, 'u_t'));
v_t = squeeze(ncread(v_t_filename, 'v_t'));

f0 = 0.909*1e-4; % Coriolis freq 
%%
%-------------------------------------------------------
%--------------------------------------------------------------------------
% line profiles at const lat
if (strcmp(config, 'nesm_2019_2020_5km'))
    lat_idx_arr = ([23]);
    lat_rho = lat_rho_nesm_5km;
    lon_rho = lon_rho_nesm_5km;
    lat_rho_vec = lat_rho_vec_nesm_5km;
    lon_rho_vec = lon_rho_vec_nesm_5km;
    depth = depth_nesm_5km;
    xyskip = 3;
else
    lat_idx_arr = ([111]);
    lat_rho = lat_rho_nesm_1km;
    lon_rho = lon_rho_nesm_1km;
    lat_rho_vec = lat_rho_vec_nesm_1km;
    lon_rho_vec = lon_rho_vec_nesm_1km;
    xyskip = 15;
    if(strcmp(config, 'nesm_2019_2020'))
        depth = depth_nesm_1km;
    elseif(strcmp(config, 'sm3_2019_2020'))
        depth = depth_sm3_1km;
    end
end
lat_arr = zeros(length(lat_idx_arr), 1);
for i = 1:length(lat_idx_arr)
        lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------

%-------------------------------------------------------
for nt = Nt-1:Nt-1
    sprintf('for loop indx = %d', nt)
    vort_t_nt = indxRange(1) - nt0 + nt;
    sprintf('vort_t arr indx = %d', vort_t_nt)
    sprintf('fig indx = %d', vort_t_nt + nt0 - 1)
        %%
        figure1 = figure(nt);
        [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        % ax1 = axes('LineWidth', 3, 'FontSize', 24);

        m_proj('miller','long', lonlim,'lat', latlim);
        % h1 = m_pcolor(lon_rho, lat_rho, log10(R));
       

        hold on;
        vort = squeeze(vort_t(vort_t_nt, :, :))./f0;
        u_hslice = squeeze(u_t(vort_t_nt, :, :));
        v_hslice = squeeze(v_t(vort_t_nt, :, :));
        zMin = -0.5;  
        zMax = 0.5; 
        h1 = m_image(lon_rho_vec, lat_rho_vec, vort);
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
       
        
        m_quiver(lon_rho(1:xyskip:end, 1:xyskip:end), lat_rho(1:xyskip:end, 1:xyskip:end), ...
         u_hslice(1:xyskip:end, 1:xyskip:end), v_hslice(1:xyskip:end, 1:xyskip:end), ...
         'color',[0 0 0]);

        hold on;
        %-----------------
        % ref arrow 
        arrow_scale_factor = 1;
        lon_pos = lonlim(2) - 0.8;
        lat_pos = latlim(1) + 0.3;
        
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
        %-----------------
        [h2, c2] = m_contourf(lon_rho, lat_rho, -depth, [vlevel vlevel], ...
            'LineWidth', 3, 'EdgeColor', [0 0 0], 'FaceColor', [1 1 1], 'FaceAlpha', 0);
        
        m_grid('tickdir','in', ...
       'xtick',([-64.99, -64 -63 -62 -61]),...  % longitude   
       'xticklabel',{'65°W', '64°W','63°W','62°W','61°W'}, ... % name longitude ticks as you want
       'ytick',([38 39]), ... % latitude        
       'yticklabel',{'38°N','39°N'}); % name latitude ticks as you want;
                
        % --- Plot horizontal section lines (constant lat) ---
        for i = 1:length(lat_arr)
            m_line([min(lon_rho_vec(:)), max(lon_rho_vec(:))], [lat_arr(i), lat_arr(i)], ...
                'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
        end

        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        title(ax1, string(date));
        set(figure1, 'Visible', 'off'); % stop pop-ups
        figname  = [plots_path, './', config, '_vort'];
        
        figname = strcat(figname, '_vlevel_', string(vlevel));
        figname = strcat(figname, '_nt_', string(indxRange(1) + nt - 1));
        % pdf
        figname = strcat(figname, '.pdf');
        exportgraphics(figure1, figname); % remove extra white space, 2022a and above

        close all;
end
% 
%%

%------------------------------------

