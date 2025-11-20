%-------------------------------------------------------
% plot averaged rho vlsice over given timeframe (typically 7 days)
% ref: Winter Mixed Layer Restratification Induced by Vertical Eddy
% Buoyancy Flux in the Labrador Sea, Li et al., GRL
% vertical eddy buoyancy flux (vebf) = <w'b'>

% Goal: to obtain b'w', with b = -g*(rho - rho0)/rho0; 
% primes denote deviations from a time average

% plot b'w'(z, t) 
% plot <b'w'>_z -> z-avgeraged time series 
%%
%-------------------------------------------------------
start_paths
%-------------------------------------------------------
indxRange = 952:3877; % What time indices do you need?
% indxRange = 2288:2416; % Oct. 15-30, 2019
% indxRange = 3024:3152; % Jan. 15-30, 2020

nt0=indxRange(1);
[~, Nt] = size(indxRange);
%-------------------------------------------------------
%-------------------------------------------------------
%%
lon_idx_arr = [150, 172];
X = (repmat(lon_rho_vec, NumLayers, 1));
%%
%-------------------------------------------------------
filename = strcat('nesm_2019_2020_KmKe_vslice_const_lat_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');

hrs_vslice_arr = (ncread(filename, 'hrs_vslice'));
vrs_vslice_arr = (ncread(filename, 'vrs_vslice'));
kmke_vslice_arr = (ncread(filename, 'kmke_vslice'));

lat_arr = (ncread(filename, 'lat'));
lat_idx_arr = (ncread(filename, 'lat_idx'));
%-------------------------------------------------------
%-------------------------------------------------------


%%
% mkdir vort_plots;
for i = 1:length(lat_idx_arr)
    lat_idx = (lat_idx_arr(i)); % find(abs(lat_rho_vec-lat_arr(1))<1e-3); % find idx of lat in lat_rho_vec
    Z = (squeeze(zr(:,lat_idx,:)));
    %-------------------------------------
    %% plotting HRS vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    hrs_vslice = squeeze(hrs_vslice_arr(:, i, :));
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, hrs_vslice); 
    colormap(cmocean('haline')); colorbar; 
    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
    
    for j = 1:length(lon_idx_arr)
            xline(squeeze(X(:, lon_idx_arr(j))));
    end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])

    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    figname  = [plots_path, 'rs_plots/hrs/nesm_2019_2020_hrs_vslice'];
    figname = strcat(figname, '_lat_', string(lat_arr(i)));
    figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
    string(indxRange(Nt)));
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %
    % save for all i
    % if(i == 1)
        % for j = 1:length(lon_idx_arr)
        %     lon_idx = lon_idx_arr(j);
        %     % save file to compare/vebf
        %     data = [squeeze(Z(:, lon_idx)), squeeze(hrs_vslice(:, lon_idx))];
        %     filename = '../../../compare/KmKe/vslice/KmKe/nesm_2019_2020_KmKe';
        %     filename = strcat(filename, '_lat_', string(lat_arr(i)));
        %     filename = strcat(filename, '_lon_', string(X(1, lon_idx)));
        %     filename = strcat(filename, '_nt_', string(indxRange(1)), '_', ...
        %     string(indxRange(Nt)));
        %     filename = strcat(filename, '.asc');
        %     dlmwrite(filename, data, 'delimiter', '\t');
        % end
    % end
    %-------------------------------------
    %% plotting VRS vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    vrs_vslice = squeeze(vrs_vslice_arr(:, i, :));
    % KmKe_vslice = hrs_vslice + vrs_vslice;
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, vrs_vslice); 
    colormap(cmocean('haline')); colorbar; 
    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
    
    for j = 1:length(lon_idx_arr)
            xline(squeeze(X(:, lon_idx_arr(j))));
    end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])

    hold on;
    
    set(figure1, 'Visible', 'off'); % stop pop-ups
    figname  = [plots_path, 'rs_plots/vrs/nesm_2019_2020_vrs_vslice'];
    figname = strcat(figname, '_lat_', string(lat_arr(i)));
    figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
    string(indxRange(Nt)));
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %-------------------------------------
    %% plotting KmKe vslice

    figure1 = figure(i);
    % [latlim, lonlim] = geoquadline(lat_rho, lon_rho);
    % Create axes
    ax1 = axes('Parent', figure1, 'YMinorTick','on',...
        'LineWidth',3,...
        'FontSize',24);
    
    hold on;
    kmke_vslice = squeeze(kmke_vslice_arr(:, i, :));
    % KmKe_vslice = hrs_vslice + vrs_vslice;
    zMin = -5e-7; 
    zMax = -zMin; % -zMin; % max(max(~isinf(w_vslice)));        % h1 = image(lon_rho_vec, z_depth_vec, w_vslice);
    pcolor(ax1, X, Z, kmke_vslice); 
    % colormap(cmocean('haline')); colorbar;
    colormap('parula'); cb = colorbar;
    cb.FontSize = 14;

    clim([zMin zMax]);
    shading interp;
    set(ax1,'Color', [1 1 1])
    hold on;
    
    % for j = 1:length(lon_idx_arr)
    %         xline(squeeze(X(:, lon_idx_arr(j))));
    % end
    
    xticks([-64.99 -64 -63 -62 -61])  % longitude   
    xticklabels({'65°W', '64°W','63°W','62°W','61°W'}) % name longitude ticks as you want
    yticks([-5000 -4000 -3000 -2000 -1000])  % height from surface in m   
    yticklabels({'5000', '4000', '3000', '2000', '1000'}) % name y ticks as you want
    % ax1.TitleHorizontalAlignment = 'left'; % left makes it come to center
    ylim([-5000 0])

    hold on;
    
    % title(ax1, string(date));

    set(figure1, 'Visible', 'off'); % stop pop-ups
    figname  = [plots_path, 'rs_plots/KmKe/nesm_2019_2020_KmKe_vslice'];
    figname = strcat(figname, '_lat_', string(lat_arr(i)));
    figname = strcat(figname, '_nt_', string(indxRange(1)), '_', ...
    string(indxRange(Nt)));
    exportgraphics(figure1, strcat(figname, '.pdf'), 'ContentType', 'vector', 'Resolution',300); % remove extra white space, 2022a and above
    % export_fig(figure1, figname, '-eps','-transparent', '-r300'); % 
    close all;
    %-------------------------------------
end
%%
