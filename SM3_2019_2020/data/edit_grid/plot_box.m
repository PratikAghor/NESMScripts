function []=plot_box(lat_vec, lon_vec, lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max, color)
    if nargin < 5
        color = 'r'; % Default color is red if not specified
    end

    m_plot([lon_vec(lon_idx_min) lon_vec(lon_idx_max)],[lat_vec(lat_idx_min) lat_vec(lat_idx_min)], ...
    color, 'LineWidth', 3);
    m_plot([lon_vec(lon_idx_min) lon_vec(lon_idx_max)],[lat_vec(lat_idx_max) lat_vec(lat_idx_max)], ...
    color, 'LineWidth', 3);

    m_plot([lon_vec(lon_idx_min) lon_vec(lon_idx_min)],[lat_vec(lat_idx_min) lat_vec(lat_idx_max)], ...
    color, 'LineWidth', 3);
    m_plot([lon_vec(lon_idx_max) lon_vec(lon_idx_max)],[lat_vec(lat_idx_min) lat_vec(lat_idx_max)], ...
    color, 'LineWidth', 3);
end