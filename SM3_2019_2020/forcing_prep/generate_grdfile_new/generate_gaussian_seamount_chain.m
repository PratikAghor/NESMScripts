function [h_new] = generate_gaussian_seamount_chain(latr, lonr, hmax, h)
% works for Ny x Nx = 278 x 300 
%------------------------
% Generate multiple Gaussian seamounts
% Creates a bathymetry with three Gaussian seamounts: Atlantis II, Kelvin, and Gosnold
%
% Inputs:
%   latr        - Latitude grid (vector)
%   lonr        - Longitude grid (vector)
%   hmax        - Maximum depth allowed
%   h           - Original bathymetry
%   nesm_grdfile - Path to NESM_grd.nc file for reference bathymetry
%
% Output:
%   h_new       - New bathymetry with three Gaussian seamounts
%
%------------------------
% Get original bathymetry size
[ny, nx] = size(h);
nesm_h = h;
% Initialize new bathymetry with original
h_new = h;
% Define seamount boxes (indices)
% Atlantis II seamount box
atl_lat_idx_min = 55;
atl_lat_idx_max = 155;
atl_lon_idx_min = 120;
atl_lon_idx_max = 195;

% Kelvin seamount idx box
kelv_lat_idx_min = 120;
kelv_lat_idx_max = 180;
kelv_lon_idx_min = 60;
kelv_lon_idx_max = 120;

% Gosnold seamount idx box
gosn_lat_idx_min = 10;
gosn_lat_idx_max = 130;
gosn_lon_idx_min = 195;
gosn_lon_idx_max = 280;

%% ------------------- KELVIN SEAMOUNT -------------------
%% ------------------- GOSNOLD SEAMOUNT -------------------
% Flatten Kelvin and Gosnold seamounts
[ny, nx] = size(h);
for j = 1:nx
    for i = 1:ny
        % Flatten Kelvin seamount
        if((j>=kelv_lon_idx_min && j<=kelv_lon_idx_max) ...
                && (i>=kelv_lat_idx_min && i <= kelv_lat_idx_max))
           h_new(i, j) = h(kelv_lat_idx_max-2, kelv_lon_idx_min+2);
        end
        
        % Flatten Gosnold seamount
        if((j>=gosn_lon_idx_min && j<=gosn_lon_idx_max) ...
                && (i>=gosn_lat_idx_min && i <= gosn_lat_idx_max))
           h_new(i, j) = h(gosn_lat_idx_max-2, gosn_lon_idx_min+2);
        end
    end
end

%% ------------------- ATLANTIS II SEAMOUNT -------------------
% Define transition zone width (buffer)
buffer = 20;

% Determine box with buffer for the Atlantis II seamount
ext_lat_min = max(1, atl_lat_idx_min - buffer);
ext_lat_max = min(ny, atl_lat_idx_max + buffer);
ext_lon_min = max(1, atl_lon_idx_min - buffer);
ext_lon_max = min(nx, atl_lon_idx_max + buffer);

boundary_depths = [];

% Sample only just outside the Atlantis II box
for j = atl_lon_idx_min-buffer:atl_lon_idx_max+buffer
    j = max(1, min(nx, j)); % column index bounds check
    row1 = max(1, min(ny, atl_lat_idx_min - buffer));
    row2 = max(1, min(ny, atl_lat_idx_max + buffer));
    boundary_depths = [boundary_depths; h_new(row1, j)];
    boundary_depths = [boundary_depths; h_new(row2, j)];
end

for i = atl_lat_idx_min-buffer:atl_lat_idx_max+buffer
    i = max(1, min(ny, i));
    col1 = max(1, min(nx, atl_lon_idx_min - buffer));
    col2 = max(1, min(nx, atl_lon_idx_max + buffer));
    boundary_depths = [boundary_depths; h_new(i, col1)];
    boundary_depths = [boundary_depths; h_new(i, col2)];
end
avg_boundary_depth = mean(boundary_depths);

% Define the Gaussian peak parameters for Atlantis II
min_depth = 1645;           % The shallowest depth at peak of seamount (summit)
atl_loncentre = -63.23;     % Longitude center
atl_latcentre = 38.5;       % Latitude center
sx = 0.03;                  % x variance
sy = 0.03;                  % y variance
base_depth = avg_boundary_depth; 
height = base_depth - min_depth; % Height from base to peak

% Pre-calculate the central region for the pure Gaussian
central_region = false(ny, nx);
for j = atl_lon_idx_min:atl_lon_idx_max
    for i = atl_lat_idx_min:atl_lat_idx_max
        central_region(i, j) = true;
    end
end

% Apply Atlantis II Gaussian seamount
h_new = apply_gaussian_seamount(h_new, h, latr, lonr, central_region, ...
                                atl_lat_idx_min, atl_lat_idx_max, ...
                                atl_lon_idx_min, atl_lon_idx_max, ...
                                buffer, atl_latcentre, atl_loncentre, ...
                                base_depth, height, sx, sy);

%% ------------------- KELVIN SEAMOUNT -------------------
% Determine boundary depths for Kelvin
kelv_boundary_depths = [];

% Sample just outside the Kelvin box
for j = kelv_lon_idx_min-buffer:kelv_lon_idx_max+buffer
    j = max(1, min(nx, j)); % column bounds
    row1 = max(1, min(ny, kelv_lat_idx_min - buffer));
    row2 = max(1, min(ny, kelv_lat_idx_max + buffer));
    kelv_boundary_depths = [kelv_boundary_depths; h_new(row1, j)];
    kelv_boundary_depths = [kelv_boundary_depths; h_new(row2, j)];
end
for i = kelv_lat_idx_min-buffer:kelv_lat_idx_max+buffer
    i = max(1, min(ny, i));
    col1 = max(1, min(nx, kelv_lon_idx_min - buffer));
    col2 = max(1, min(nx, kelv_lon_idx_max + buffer));
    kelv_boundary_depths = [kelv_boundary_depths; h_new(i, col1)];
    kelv_boundary_depths = [kelv_boundary_depths; h_new(i, col2)];
end
kelv_avg_boundary_depth = mean(kelv_boundary_depths);

% Define the Gaussian peak parameters for Kelvin seamount
% extract parameters from reference bathymetry
kelv_region = nesm_h(kelv_lat_idx_min:kelv_lat_idx_max, kelv_lon_idx_min:kelv_lon_idx_max);
kelv_min_depth = min(min(kelv_region));
[kelv_lat_idx, kelv_lon_idx] = find(kelv_region == kelv_min_depth, 1, 'first');
kelv_lat_idx = kelv_lat_idx + kelv_lat_idx_min - 1;
kelv_lon_idx = kelv_lon_idx + kelv_lon_idx_min - 1;
kelv_latcentre = latr(kelv_lat_idx);
kelv_loncentre = lonr(kelv_lon_idx);

kelv_sx = 0.02;  % x variance - slightly narrower than Atlantis II
kelv_sy = 0.02;  % y variance
kelv_base_depth = kelv_avg_boundary_depth;
kelv_height = kelv_base_depth - kelv_min_depth;

% Pre-calculate the central region for Kelvin
kelv_central_region = false(ny, nx);
for j = kelv_lon_idx_min:kelv_lon_idx_max
    for i = kelv_lat_idx_min:kelv_lat_idx_max
        kelv_central_region(i, j) = true;
    end
end

% Apply Kelvin Gaussian seamount
h_new = apply_gaussian_seamount(h_new, h, latr, lonr, kelv_central_region, ...
                                kelv_lat_idx_min, kelv_lat_idx_max, ...
                                kelv_lon_idx_min, kelv_lon_idx_max, ...
                                buffer, kelv_latcentre, kelv_loncentre, ...
                                kelv_base_depth, kelv_height, kelv_sx, kelv_sy);

%% ------------------- GOSNOLD SEAMOUNT -------------------
buffer = 15;
% Determine boundary depths for Gosnold
gosn_boundary_depths = [];

% Sample just outside the Gosnold box
for j = gosn_lon_idx_min-buffer:gosn_lon_idx_max+buffer
    j = max(1, min(nx, j)); % column bounds
    row1 = max(1, min(ny, gosn_lat_idx_min - buffer));
    row2 = max(1, min(ny, gosn_lat_idx_max + buffer));
    gosn_boundary_depths = [gosn_boundary_depths; h_new(row1, j)];
    gosn_boundary_depths = [gosn_boundary_depths; h_new(row2, j)];
end
for i = gosn_lat_idx_min-buffer:gosn_lat_idx_max+buffer
    i = max(1, min(ny, i));
    col1 = max(1, min(nx, gosn_lon_idx_min - buffer));
    col2 = max(1, min(nx, gosn_lon_idx_max + buffer));
    gosn_boundary_depths = [gosn_boundary_depths; h_new(i, col1)];
    gosn_boundary_depths = [gosn_boundary_depths; h_new(i, col2)];
end

gosn_avg_boundary_depth = mean(gosn_boundary_depths);

% Define the Gaussian peak parameters for Gosnold seamount
% extract parameters from reference bathymetry
gosn_region = nesm_h(gosn_lat_idx_min:gosn_lat_idx_max, gosn_lon_idx_min:gosn_lon_idx_max);
gosn_min_depth = min(min(gosn_region));
[gosn_lat_idx, gosn_lon_idx] = find(gosn_region == gosn_min_depth, 1, 'first');
gosn_lat_idx = gosn_lat_idx + gosn_lat_idx_min - 1;
gosn_lon_idx = gosn_lon_idx + gosn_lon_idx_min - 1;
gosn_latcentre = latr(gosn_lat_idx);
gosn_loncentre = lonr(gosn_lon_idx);

gosn_sx = 0.03;  % x variance
gosn_sy = 0.03;  % y variance
gosn_base_depth = gosn_avg_boundary_depth;
gosn_height = gosn_base_depth - gosn_min_depth;

% Pre-calculate the central region for Gosnold
gosn_central_region = false(ny, nx);
for j = gosn_lon_idx_min:gosn_lon_idx_max
    for i = gosn_lat_idx_min:gosn_lat_idx_max
        gosn_central_region(i, j) = true;
    end
end

% Apply Gosnold Gaussian seamount
h_new = apply_gaussian_seamount(h_new, h, latr, lonr, gosn_central_region, ...
                                gosn_lat_idx_min, gosn_lat_idx_max, ...
                                gosn_lon_idx_min, gosn_lon_idx_max, ...
                                buffer, gosn_latcentre, gosn_loncentre, ...
                                gosn_base_depth, gosn_height, gosn_sx, gosn_sy);


%---------------------------------------------------------------------
%% Reapply Clean Atlantis II Gaussian
% Recompute the Gaussian height field
% reset buffer to -5
buffer =-5;
% Determine box with buffer for the Atlantis II seamount
ext_lat_min = max(1, atl_lat_idx_min - buffer);
ext_lat_max = min(ny, atl_lat_idx_max + buffer);
ext_lon_min = max(1, atl_lon_idx_min - buffer);
ext_lon_max = min(nx, atl_lon_idx_max + buffer);

% Extract surrounding bathymetry at the boundary for blending reference
boundary_depths = [];

% Sample only just outside the Atlantis II box
for j = atl_lon_idx_min-buffer:atl_lon_idx_max+buffer
    j = max(1, min(nx, j)); % column index bounds check
    row1 = max(1, min(ny, atl_lat_idx_min - buffer));
    row2 = max(1, min(ny, atl_lat_idx_max + buffer));
    boundary_depths = [boundary_depths; h_new(row1, j)];
    boundary_depths = [boundary_depths; h_new(row2, j)];
end

for i = atl_lat_idx_min-buffer:atl_lat_idx_max+buffer
    i = max(1, min(ny, i)); % row index bounds check
    col1 = max(1, min(nx, atl_lon_idx_min - buffer));
    col2 = max(1, min(nx, atl_lon_idx_max + buffer));
    boundary_depths = [boundary_depths; h_new(i, col1)];
    boundary_depths = [boundary_depths; h_new(i, col2)];
end
avg_boundary_depth = mean(boundary_depths);

% Define the Gaussian peak parameters for Atlantis II
min_depth = 1645;           % The shallowest depth at peak of seamount (summit)
atl_loncentre = -63.23;     % Longitude center
atl_latcentre = 38.5;       % Latitude center
sx = 0.03;                  % x variance
sy = 0.03;                  % y variance
base_depth = avg_boundary_depth; 
height = base_depth - min_depth; % Height from base to peak

% Pre-calculate the central region for the pure Gaussian
central_region = false(ny, nx);
for j = atl_lon_idx_min:atl_lon_idx_max
    for i = atl_lat_idx_min:atl_lat_idx_max
        central_region(i, j) = true;
    end
end

% Apply Atlantis II Gaussian seamount
h_new = apply_gaussian_seamount(h_new, h, latr, lonr, central_region, ...
                                atl_lat_idx_min, atl_lat_idx_max, ...
                                atl_lon_idx_min, atl_lon_idx_max, ...
                                buffer, atl_latcentre, atl_loncentre, ...
                                base_depth, height, sx, sy);
%--------------------------------------------------------------------------
% Restore original bathymetry outside all extended regions with buffers
buffer = 0;

ext_atl_lat_min = max(1, atl_lat_idx_min - buffer);
ext_atl_lat_max = min(ny, atl_lat_idx_max + buffer);
ext_atl_lon_min = max(1, atl_lon_idx_min - buffer);
ext_atl_lon_max = min(nx, atl_lon_idx_max + buffer);

ext_kelv_lat_min = max(1, kelv_lat_idx_min - buffer);
ext_kelv_lat_max = min(ny, kelv_lat_idx_max + buffer);
ext_kelv_lon_min = max(1, kelv_lon_idx_min - buffer);
ext_kelv_lon_max = min(nx, kelv_lon_idx_max + buffer);

ext_gosn_lat_min = max(1, gosn_lat_idx_min - buffer);
ext_gosn_lat_max = min(ny, gosn_lat_idx_max + buffer);
ext_gosn_lon_min = max(1, gosn_lon_idx_min - buffer);
ext_gosn_lon_max = min(nx, gosn_lon_idx_max + buffer);

% Create a mask of all modified areas
modified_mask = false(ny, nx);
modified_mask(ext_atl_lat_min:ext_atl_lat_max, ext_atl_lon_min:ext_atl_lon_max) = true;
modified_mask(ext_kelv_lat_min:ext_kelv_lat_max, ext_kelv_lon_min:ext_kelv_lon_max) = true;
modified_mask(ext_gosn_lat_min:ext_gosn_lat_max, ext_gosn_lon_min:ext_gosn_lon_max) = true;

% Restore original bathymetry outside the modified regions
for j = 1:nx
    for i = 1:ny
        if ~modified_mask(i, j)
            h_new(i, j) = h(i, j);
        end
    end
end
%--------------------------------------------------------------------------
end