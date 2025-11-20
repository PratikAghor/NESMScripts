function h_out = apply_gaussian_seamount(h_in, h_orig, latr, lonr, central_region, ...
                                        lat_idx_min, lat_idx_max, ...
                                        lon_idx_min, lon_idx_max, ...
                                        buffer, latcentre, loncentre, ...
                                        base_depth, height, sx, sy)
    % Apply a Gaussian seamount to the bathymetry
    %
    % Inputs:
    %   h_in          - Current bathymetry
    %   h_orig        - Original bathymetry
    %   latr, lonr    - Lat/lon grid vectors
    %   central_region - Boolean mask for central region
    %   lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max - Region indices
    %   buffer        - Buffer width for smooth transition
    %   latcentre, loncentre - Center coordinates
    %   base_depth    - Depth at base of seamount
    %   height        - Height of seamount from base to peak
    %   sx, sy        - Variance parameters
    %
    % Output:
    %   h_out         - Updated bathymetry with Gaussian seamount
%--------------------------------------------------------------------    
    h_out = h_in;
    [ny, nx] = size(h_in);
    
    % Determine extended region with buffer
    ext_lat_min = max(1, lat_idx_min - buffer);
    ext_lat_max = min(ny, lat_idx_max + buffer);
    ext_lon_min = max(1, lon_idx_min - buffer);
    ext_lon_max = min(nx, lon_idx_max + buffer);
    
    % For each point in the extended region
    for j = ext_lon_min:ext_lon_max
        for i = ext_lat_min:ext_lat_max
            local_lat = latr(i);
            local_lon = lonr(j);
            
            % Calculate pure Gaussian function
            gaussian_val = exp(-(local_lon - loncentre)^2/(sx) - (local_lat - latcentre)^2/(sy));
            
            dist_from_edge = 0;
            
            % Check if we are in the central region or the buffer zone
            if central_region(i, j)
                % In central region - full Gaussian effect
                dist_from_edge = 0;
            else
                % In buffer zone - calculate distance to central region
                if i < lat_idx_min
                    dist_lat = lat_idx_min - i;
                elseif i > lat_idx_max
                    dist_lat = i - lat_idx_max;
                else
                    dist_lat = 0;
                end
                
                if j < lon_idx_min
                    dist_lon = lon_idx_min - j;
                elseif j > lon_idx_max
                    dist_lon = j - lon_idx_max;
                else
                    dist_lon = 0;
                end
                
                dist_from_edge = max(dist_lat, dist_lon);
            end
            
            % Apply smooth transition at the boundary
            if dist_from_edge == 0
                % We're inside the central region - apply full Gaussian
                gaussian_depth = base_depth - height * gaussian_val;
                h_out(i, j) = gaussian_depth;
            elseif dist_from_edge <= buffer
                % We're in the transition zone - use a smooth cosine blend
                blend_factor = 0.5 * (1 + cos(pi * dist_from_edge / buffer));
                
                % Calculate the Gaussian depth
                gaussian_depth = base_depth - height * gaussian_val;
                
                % Use original bathymetry for blending to ensure perfect transition
                original_depth = h_orig(i, j);
                
                % Blend between Gaussian depth and original depth
                h_out(i, j) = blend_factor * gaussian_depth + (1 - blend_factor) * original_depth;
            end
            % Outside buffer zone: h_out already contains original bathymetry
        end
    end
end