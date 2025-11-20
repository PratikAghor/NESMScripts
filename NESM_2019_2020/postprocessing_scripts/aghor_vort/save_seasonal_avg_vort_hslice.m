% --------------------------------------------------------
% Save seasonal mean vorticity horizontal slices (depth = vlevel)
% Output shape: (4, Ny, Nx)
% --------------------------------------------------------
clear all; clc;
start_paths;

% Parameters
vlevel = -500;    % depth level [m]
coef = 1;          % scaling factor if needed
seasons = {'JJA', 'SON', 'DJF', 'MA'};
season_inds = {
    1200:1935,     % JJA
    1936:2663,     % SON
    2664:3391,     % DJF
    3392:3877      % MA
};

% Initialize
% sample_file = sprintf('%sNESM_2019_2020_avg.%05d.nc', dirHR, season_inds{1}(1));
%[~, ~, ~, sample_vort] = get_vort(sample_file, gridfile, season_inds{1}(1), vlevel, coef);
%[Ny, Nx] = size(sample_vort);

vort_seasonal_avg = zeros(4, Ny, Nx);  % shape = (season, y, x)

% Loop through each season
for s = 1:4
    inds = season_inds{s};
    Nt = length(inds);
    vort_sum = zeros(Ny, Nx);

    for nt = 1:Nt
        n = inds(nt);
        fname = sprintf('%sNESM_2019_2020_avg.%05d.nc', dirHR, n);
        fprintf('Loading %s (%s, index %d of %d)\n', fname, seasons{s}, nt, Nt);
        [~, ~, ~, vort_slice] = get_vort(fname, gridfile, tindex, vlevel, coef);
        vort_sum = vort_sum + vort_slice;
    end

    vort_seasonal_avg(s, :, :) = vort_sum / Nt;
end

% --------------------------------------------------------
% Write output to NetCDF
% --------------------------------------------------------
outfile = sprintf('nesm_2019_2020_seasonal_avg_vort_vlevel_%d.nc', (vlevel));
ncid = netcdf.create(outfile, 'CLOBBER');

% Define dimensions
dim_x = netcdf.defDim(ncid, 'Nx', Nx);
dim_y = netcdf.defDim(ncid, 'Ny', Ny);
dim_season = netcdf.defDim(ncid, 'season', 4);

% Define variable (dimensions: [Nx, Ny, season])
vort_varid = netcdf.defVar(ncid, 'vort_t_avg', 'double', [dim_season, dim_y, dim_x]);
netcdf.putAtt(ncid, vort_varid, 'description', ...
    sprintf('Seasonal average vorticity slices at depth = %dm', vlevel));
netcdf.putAtt(ncid, vort_varid, 'units', 's^-1');

% (Optional) Add season name variable
dim_strlen = netcdf.defDim(ncid, 'strlen', 3);
season_str_id = netcdf.defVar(ncid, 'season_names', 'char', [dim_strlen, dim_season]);
netcdf.putAtt(ncid, season_str_id, 'description', 'Season abbreviations');

% End define mode
netcdf.endDef(ncid);

% Write data (permute for NetCDF: [season, Ny, Nx])
netcdf.putVar(ncid, vort_varid, vort_seasonal_avg);
netcdf.putVar(ncid, season_str_id, ['JJA'; 'SON'; 'DJF'; 'MA_']');

netcdf.close(ncid);

disp('Done! Saved seasonal-averaged vorticity to NetCDF.');

