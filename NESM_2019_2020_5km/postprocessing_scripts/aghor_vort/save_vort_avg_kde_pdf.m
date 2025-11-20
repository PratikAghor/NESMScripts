%---------------------------------------------------
% 1. calculate vort data on zlevels
% 2. z-Average vort/f0 from -1850 m to bottom
% 3. calculate kde pdf of the resulting dimensionless vort ([Ny x Nx])
% 4. Average these pdfs over the specified timespan
% 5. Save the mean vort_bins and pdf in a netcdf file
% Author: Pratik Aghor
%---------------------------------------------------
clc;
clear all;
close all;
start_paths;  % make sure this sets your dirHR, gridfile, coef, Nx, Ny, N, etc.

%-------------------- CONFIG ------------------------%
% indxRange = 952:3877;  % Time indices for the loop
indxRange=3024:3152; % Jan 15 - Jan 30, 2020

f0 = 0.909e-4;         % Coriolis parameter (s^-1)
z_bottom = -5500;      % Bottom depth (m)
z_threshold = -1850;   % Depth to start vertical averaging (m)
vort_bins = linspace(-1, 1, 400);  % KDE evaluation points (dimensionless vort range)

%-------------------- ZLEVS SETUP --------------------%
zlevs = linspace(0, z_bottom, N); % vertical depths [0 (surface) to -5500m]
vlevel_arr = zlevs;

% Indices of levels deeper than or equal to -1850 m
deep_inds = find(vlevel_arr <= z_threshold);

Nt = length(indxRange);

pdf_matrix = zeros(Nt, length(vort_bins));  % store KDE for each timestep
%--------------------------------
for nt = 1:Nt
    vort_nt = indxRange(nt);
    fprintf('Processing time index %d (%d of %d)\n', vort_nt, nt, Nt);

    fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], vort_nt);
    hisfile = fname;

    vort_t_slices = zeros(N, Ny, Nx);

    %--- Compute full vertical slice of vorticity ---%
    for k = 1:N
        vlevel = vlevel_arr(k);
        [~, ~, ~, vort_t_slices(k, :, :)] = get_vort(hisfile, gridfile, 1, vlevel, coef);
    end

    %--- Average over depth: z = -1850 m to bottom ---%
    vort_deep = vort_t_slices(deep_inds, :, :);
    vortf_deep = vort_deep ./ f0;  % dimensionless vort

    % Vertical mean over selected depth range
    vortf_vertmean = squeeze(mean(vortf_deep, 1, 'omitnan'));  % [Ny x Nx];  % [Ny x Nx]

    % Flatten field to vector & remove NaNs
    vortf_vec = vortf_vertmean(:);
    vortf_vec = vortf_vec(~isnan(vortf_vec));

    %--- Kernel density estimate ---%
    [pdf_vals, ~] = ksdensity(vortf_vec, vort_bins);
    pdf_matrix(nt, :) = pdf_vals;
end
%------------------------------------
%--- Average KDE over time ---%
mean_pdf = mean(pdf_matrix, 1);

%--- Save results to NetCDF ---%
nc_filename = 'nesm_2019_2020_5km_dimless_vort_kde_avg_1850m_to_bottom.nc';

% Create variables and write data
nccreate(nc_filename, 'vort_bins', 'Dimensions', {'vort_bin_dim', length(vort_bins)});
nccreate(nc_filename, 'mean_pdf', 'Dimensions', {'vort_bin_dim', length(mean_pdf)});

ncwrite(nc_filename, 'vort_bins', vort_bins);
ncwrite(nc_filename, 'mean_pdf', mean_pdf);

% Add global attributes
ncwriteatt(nc_filename, '/', 'title', 'Time-averaged KDE PDF of dimensionless vorticity');
ncwriteatt(nc_filename, '/', 'depth_range', sprintf('from %.0f m to bottom (%.0f m)', z_threshold, z_bottom));
ncwriteatt(nc_filename, '/', 'units_vort_bins', 'dimensionless vorticity (\zeta/f0)');
ncwriteatt(nc_filename, '/', 'description', ...
    'Kernel Density Estimate (KDE) of spatially averaged vorticity normalized by f0, averaged over time');

fprintf('Saved KDE PDF to NetCDF file: %s\n', nc_filename);
%------------------------------------
%--- Plot final averaged PDF ---%
% figure;
% plot(vort_bins, mean_pdf, 'k', 'LineWidth', 2);
% xlabel('\zeta/f_0 (dimensionless)', 'FontSize', 14);
% ylabel('Mean PDF', 'FontSize', 14);
% title('Time-averaged PDF of \zeta/f_0 (1850 m to bottom)', 'FontSize', 16);
% grid on;
%------------------------------------
