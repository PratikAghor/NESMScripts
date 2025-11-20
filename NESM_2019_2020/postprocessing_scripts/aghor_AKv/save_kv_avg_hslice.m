%
% save time-avg kv hslice at specific depth
% and time-stamps
%
%-------------------------------------------------------
start_paths
type='r'; % see get_hslice
% tindex=1; 
%-------------------------------------------------------
indxRange = 952:3877; 
% indxRange = 3024:3151; % Jan 15 - Jan 30, 2020
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
t_arr = zeros(Nt, 1);

%% plot options
vlevel_skip = 25;
legendCell = [];
vlevel = -4000; % vlevel < 0 => depth;
Nz = 1; % length(vlevel_arr);
kv_avg_hslice = zeros(Ny, Nx);
%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    t_arr(nt) = time;
    
    AKt = pagetranspose(ncread(fname, 'AKt')); % size [Nx, Ny, Nz+1]
    AKt = shiftdim(AKt, 2); % now it is (Nz+1) x Ny x Nx, like I prefer, zw also has (Nz+1) x Ny x Nx size

    % Interpolate AKt to rho points
    AKt_rho = 0.5 * (AKt(1:N,:,:) + AKt(2:N+1,:,:));

    if(nt == 1)
                size(AKt_rho)
    end

    kv_hslice = vinterp(AKt_rho, zr, vlevel);
    kv_avg_hslice = kv_avg_hslice + kv_hslice; % get hslice at vlevel
end

kv_avg_hslice = kv_avg_hslice./Nt; % time avg
%%
%----------------------------------------------------------------------------
% save
% save data into a netcdf file
filename = strcat('nesm_2019_2020_kv_avg_hslice_vlevel_', string(vlevel) , '_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)),'.nc');

ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

% time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
% netcdf.putAtt(ncid, time_varid, 'description', 'time array');
% netcdf.putAtt(ncid, time_varid, 'units', 's');
% netcdf.putAtt(ncid, time_varid, 'array dimensions', size(t_arr));

kv_varid = netcdf.defVar(ncid, 'kv_avg_hslice', 'double', [y_len x_len]);
netcdf.putAtt(ncid, kv_varid, 'description', 't-avg of kv hslice at const vlevel (Ny x Nx)');
netcdf.putAtt(ncid, kv_varid, 'units', 'm^2/s}');
netcdf.putAtt(ncid, kv_varid, 'array dimensions', size(kv_avg_hslice));
%close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, time_varid, t_arr); % start from 0
netcdf.putVar(ncid, kv_varid, kv_avg_hslice);
% close netcdf file
netcdf.close(ncid);
disp('Done creating kv_avg_hslice at const vlevel!')
%--------------------------------------------------------
%%

   
