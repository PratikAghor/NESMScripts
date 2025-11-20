clear all
addpath('../../aghor_extras/m_map/'); 
plots_path = './';

cases = {
    {'../fig_3/NESM_grd.nc', 'AtlantisII_1km', 'nesm_bathy.pdf'};
    {'../fig_3/NESM_grd_5km.nc', 'AtlantisII_5km', 'nesm_bathy_5km.pdf'};
    {'../fig_3/SM3_grd.nc', 'AtlantisII_1km', 'sm3_bathy.pdf'};
};

for i = 1:length(cases)
    ncfile = cases{i}{1};
    boxname = cases{i}{2};
    figname = cases{i}{3};

    lon = pagetranspose(ncread(ncfile, 'lon_rho'));
    lat = pagetranspose(ncread(ncfile, 'lat_rho'));
    depth = pagetranspose(ncread(ncfile, 'h'));

    lon_rho_vec = squeeze(lon(1, :));
    lat_rho_vec = squeeze(lat(:, 1));

    [lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(boxname);

    figure1 = figure();
    [latlim, lonlim] = geoquadline(lat, lon);
    m_proj('miller', 'long', lonlim, 'lat', latlim);
    m_contourf(lon, lat, depth);
    clim([0 5500]);

    m_grid('tickdir', 'in');
    hold on;

    plot_box(lat_rho_vec, lon_rho_vec, lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max, 'r');

    if contains(figname, 'sm1_bathy')
        h = colorbar;
        h.Label.String = '[m]';
        h.Label.FontSize = 24;
        h.FontSize = 12;
        h.YTick = [1000 2000 3000 4000 5000 5500];
        h.YTickLabel = {'1000', '2000', '3000', '4000', '5000', '5500'};
    end

    set(figure1, 'Visible', 'off');
    exportgraphics(figure1, figname, 'ContentType', 'vector');
    close(figure1);
end
