%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
indxRange=3112:3119; % Jan 26, 2020
% config = 'nesm_2019_2020';
% config = 'nesm_2019_2020_5km';
config = 'sm3_2019_2020';

nt0=indxRange(1);
% Aghor's script
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

Omega = 7.2921e-5; % rotation rate of Earth s^-1
%-------------------------------------------------------
%-------------------------------------------------------
%%
vort_t_filename = strcat(config, '_vort_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
t_arr = ncread(vort_t_filename, 'time');
lat_idx_arr = ncread(vort_t_filename, 'lat_idx_arr');
lat_arr = ncread(vort_t_filename, 'lat_arr');
vort_t_vslice = ncread(vort_t_filename, 'vort_t_vslice_const_lon'); % saved it as const_lon by mistake!
%%
%-------------------------------------------------------
if (strcmp(config, 'nesm_2019_2020_5km'))
    lat_idx_arr = ([23]);
    lat_rho = lat_rho_nesm_5km;
    lon_rho = lon_rho_nesm_5km;
    lat_rho_vec = lat_rho_vec_nesm_5km;
    lon_rho_vec = lon_rho_vec_nesm_5km;
    depth = depth_nesm_5km;
    Z = zr_const_lat_5km;
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
    Z = zr_const_lat_1km;
end

% mkdir vort_plots;
%% plotting
for nt = Nt-1:Nt-1
    sprintf('for loop indx = %d', nt)
    vort_t_nt = indxRange(1) - nt0 + nt;
    sprintf('vort_t arr indx = %d', vort_t_nt)
    sprintf('fig indx = %d', vort_t_nt + nt0 - 1)
    for i = 1:1 % length(lat_arr)
        f0 = 2*Omega*sin(lat_arr(i)); % Coriolis freq.
        
        lat_idx = lat_idx_arr(i);
        X = (repmat(lon_rho_vec, N, 1));
        %%% plotting
        figure1 = figure(nt);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        vort_vslice = squeeze(vort_t_vslice(vort_t_nt, :, i, :));
        zMin = -0.5; % min(min(vort_vslice)); 
        zMax = -zMin; % max(max(~isinf(w_vslice)));
        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
        pcolor(ax1, X, Z, (vort_vslice)./f0); 
        shading interp;
        set(ax1,'Color', [1 1 1])
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        colorbar; 
       
        % bathy
        plot(ax1, lon_rho_vec, squeeze(Z(1, :)), 'k-', 'LineWidth', 3);
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        % title(ax1, string(date));

        ax1.TitleHorizontalAlignment = 'center'; % 

        xticks([-64.99 -64 -63 -62]);
        xticklabels({'65°W','64°W','63°W','62°W'});
        yticks([-5000 -4000 -3000 -2000 -1000 0]);
        yticklabels({'5','4','3','2','1','0'});
        ylabel('Depth($\times 10^3$ m)', 'Interpreter', 'latex'); 
        ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
        ylim(ax1, [-5500 0])
    %%%
        set(figure1, 'Visible', 'off'); % stop pop-ups
        figname  = [plots_path, './', config, '_vort_vslice'];
        
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(1) + nt - 1));
        % vort_contour = strcat(vort_contour, '.pdf');
        exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        % 
        close all;
    end
end
% 
%%

%------------------------------------

