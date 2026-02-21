%---------------------------------
% save model kv at const lon 
% Author: Pratik Aghor
%---------------------------------
clear all; clc; start_paths;

%---------------------------------------------------
% indxRange = 952:3877; % entire year
% indxRange=2176:2296; % Oct 1 - Oct 15, 2019
% indxRange=2288:2416; % Oct 15 - Oct 30, 2019
indxRange = 3088:3095; % Jan 23, 2020
% indxRange=3024:3152; % Jan 15 - Jan 30, 2020
nt0=indxRange(1);
[~, Nt] = size(indxRange);

Mu = M; Lu = L-1; Mv = M-1; Lv = L;
%------------------------------------------------------
% for vslice at const_lat
% lat_idx_arr = ([111]);

% lat_arr = zeros(length(lat_idx_arr), 1); % const lat vals to save vslice at
% for i = 1:length(lat_idx_arr)
%         lat_arr(i) = lat_rho_vec(lat_idx_arr(i), 1);
% end
% lonsec = lon_rho_vec;
% X = (repmat(lon_rho_vec, N, 1));
%------------------------------------------------------
% for vslice at const lon
lon_idx_arr=([159]); % 63.15 W
% lon_idx_arr = ([171]); % 63 W
lon_arr = zeros(length(lon_idx_arr), 1);
for j = 1:length(lon_idx_arr)
        lon_arr(j) = lon_rho_vec(1, lon_idx_arr(j));
end

lonval = lon_arr(1)

latsec = lat_rho_vec;
%---------------------------------------------------
kv_t_vslice_const_lon = zeros(Nt, N, Ny, length(lon_idx_arr));
t_arr = zeros(Nt, 1);
%---------------------------------------------------
%-------------------------------------------------------
for nt = 1:Nt
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(nt));
    fprintf('Loading file: %s\n', fname);

    t_arr(nt) = ncread(fname, 'time');

    AKt = pagetranspose(ncread(fname, 'AKt')); % size [Nx, Ny, Nz+1]
    AKt = shiftdim(AKt, 2); % now it is (Nz+1) x Ny x Nx, like I prefer, zw also has (Nz+1) x Ny x Nx size

    % Interpolate AKt to rho points
    AKt_rho = 0.5 * (AKt(1:N,:,:) + AKt(2:N+1,:,:));

    if(nt == 1)
                size(AKt_rho)
    end

    for j = 1:length(lon_idx_arr)
        lon_idx = lon_idx_arr(j);
        kv_t_vslice_const_lon(nt, :, :, j) = squeeze(AKt_rho(:, :, lon_idx));
    end

    squeeze(kv_t_vslice_const_lon(nt, :, 133, 1)) % at 38.7 N, 63 W

end

%--------------------------------------
% save data into a netcdf file

% filename = strcat('nesm_2019_2020_model_kv_t_vslice_const_lon_', lonval,'_nt_', string(indxRange(1)), '_', string(indxRange(Nt)), '.nc');

filename = sprintf('nesm_2019_2020_model_kv_t_vslice_const_lon_%.2f_nt_%d_%d.nc', lonval, indxRange(1), indxRange(Nt));
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', length(lon_idx_arr));
y_len = netcdf.defDim(ncid, 'Ny', Ny);
z_len = netcdf.defDim(ncid, 'Nz', NumLayers);
t_len = netcdf.defDim(ncid, 'Nt', Nt);

netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

t_varid = netcdf.defVar(ncid, 'time', 'double', [t_len]);
netcdf.putAtt(ncid, t_varid, 'description', 'time');
netcdf.putAtt(ncid, t_varid, 'units', 's');
netcdf.putAtt(ncid, t_varid, 'array dimensions', size(t_arr));

lon_idx_varid = netcdf.defVar(ncid, 'lon_idx_arr', 'double', [x_len]);
netcdf.putAtt(ncid, lon_idx_varid, 'description', 'indices in lon_rho_vec for const lon values where vertical slice is taken');
netcdf.putAtt(ncid, lon_idx_varid, 'units', '--');
netcdf.putAtt(ncid, lon_idx_varid, 'array dimensions', size(lon_idx_arr));

lon_varid = netcdf.defVar(ncid, 'lon_arr', 'double', [x_len]);
netcdf.putAtt(ncid, lon_varid, 'description', 'const lon values where vertical slice is taken');
netcdf.putAtt(ncid, lon_varid, 'units', '--');
netcdf.putAtt(ncid, lon_varid, 'array dimensions', size(lon_arr));

w_varid = netcdf.defVar(ncid, 'kv_t_vslice_const_lon', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'timeseries of kv_t at const_lon');
netcdf.putAtt(ncid, w_varid, 'units', 'm^2/s');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(kv_t_vslice_const_lon));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, lon_idx_varid, lon_idx_arr);
netcdf.putVar(ncid, lon_varid, lon_arr);
netcdf.putVar(ncid, w_varid, kv_t_vslice_const_lon);
% close netcdf file
netcdf.close(ncid);
%----------------------------------------
% save data into a netcdf file
disp('Done saving kv_t_vslice_const_lon!')
%--------------------------------------------------------
