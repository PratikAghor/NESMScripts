%-------------------------------------------------------
%%
start_paths
%-------------------------------------------------------
Omega = 7.2921e-5; % rotation rate of Earth s^-1

indxRange = 952:3877;
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
% indxRange=3024:3152; % Jan 15 - Jan 30, 2020
% indxRange=3264:3392; % Feb 14 - Feb 29, 2020

nt0=indxRange(1);
% Aghor's script
% get avg ke as a function of height and temperature
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);

tskip = 1;
%-------------------------------------------------------
%-------------------------------------------------------
%%
vort_t_filename = strcat('nesm_2019_2020_5km_vort_t_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');
% vort_t_filename = 'nesm_2019_vort_t_1000_2902.nc';
% nt0 = 1000; % 0 if 0_999, 1000 if 1000_2902
% nt0=indxRange(1);
t_arr = ncread(vort_t_filename, 'time');
lat_idx_arr = ncread(vort_t_filename, 'lat_idx_arr');
lat_arr = ncread(vort_t_filename, 'lat_arr');
vort_t_vslice = ncread(vort_t_filename, 'vort_t_vslice_const_lon');
%%
%-------------------------------------------------------
for nt = 2136:tskip:2139
    sprintf('for loop indx = %d', nt)
    vort_t_nt = indxRange(1) - nt0 + nt;
    sprintf('vort_t arr indx = %d', vort_t_nt)
    sprintf('fig indx = %d', vort_t_nt + nt0 - 1)
    for i = 1:1 % length(lat_arr)
        f0 = 2*Omega*sin(lat_arr(i)); % Coriolis freq.
        
        lat_idx = lat_idx_arr(i);
        z_depth_vec = squeeze(zr(:, lat_idx, :)); % get depth vals at that lat
        Z = ((squeeze(zr(:,lat_idx,:))));
        X = (repmat(lon_rho_vec, NumLayers, 1));
        %%% plotting
        figure1 = figure(nt);
        % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
        % Create axes
        ax1 = axes('Parent', figure1, 'YMinorTick','on',...
            'LineWidth',3,...
            'FontSize',24);
        hold on;
        vort_vslice = squeeze(vort_t_vslice(vort_t_nt, :, i, :));
        zMin = -0.5; 
        zMax = -zMin; 
        pcolor(ax1, X, Z, (vort_vslice)./f0); 
        shading interp;
        set(ax1,'Color', [1 1 1])
        colormap(cmocean('balance')); colorbar; clim([zMin zMax]);
        colorbar; 
        date = char(string(datetime(2017, 12, 31, 23, 29, t_arr(vort_t_nt))));
        date = date(1:11);
        % title(ax1, string(date));

        ax1.TitleHorizontalAlignment = 'center'; 

        xticks([-64.99 -64 -63 -62 -61])  
        xticklabels({'65°W', '64°W','63°W','62°W','61°W'})
        yticks([-5000 -4000 -3000 -2000 -1000])  
        yticklabels({'5000', '4000', '3000', '2000', '1000'}) 
        ax1.TitleHorizontalAlignment = 'left';
        ylim([-5500 0])
        set(figure1, 'Visible', 'off'); % stop pop-ups
        figname  = [plots_path, 'vort_plots/vslice/const_lat/nesm_2019_2020_5km_vort_vslice'];
        
        figname = strcat(figname, '_lat_', string(lat_arr(i)));
        figname = strcat(figname, '_nt_', string(indxRange(1) + nt - 1));
        exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector'); % remove extra white space, 2022a and above
        close all;
    end
end
% 
%%

%------------------------------------

