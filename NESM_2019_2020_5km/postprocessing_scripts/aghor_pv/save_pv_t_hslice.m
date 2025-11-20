%
% save horizontal slices of w for selected vertical levles (vlevels)
% and time-stamps
%
%-------------------------------------------------------
start_paths
type='r'; % see get_hslice
% tindex=1; 
%-------------------------------------------------------
indxRange = 0:3877; % What time indices do you need?
init = 0; %Index of first time period
% Aghor's script
[~, Nt] = size(indxRange);
hke_arr = zeros(Nt, NumLayers); % get avg w as a function of height and time
t_arr = zeros(Nt, 1);

legendCell = [];
% vlevel_arr = ([100, 52, 39, 15, 1]);% NumLayers: -vlevel_skip: 1
vlevel = -4000; % vlevel < 0 => depth;
Nz = 1; % length(vlevel_arr);
epv_t_hslice = zeros(Nt, Ny, Nx);
%%
for nt = 1:Nt
    sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
    fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    time = ncread(fname, 'time');
    t_arr(nt) = time;
    epv_3d = ertel_aghor(hisfile,gridfile,'rho',tindex); 
    epv_t_hslice(nt, :, :) = vinterp(epv_3d, zr, vlevel); % get hslice at vlevel
end

epv_t_hslice = -g.*epv_t_hslice;

%%
%----------------------------------------------------------------------------
% save
% save data into a netcdf file
filename = strcat('nesm_2019_2020_5km_epv_t_hslice_vlevel_', string(vlevel) , '_nt_',string(indxRange(1)), '_', ...
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

time_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, time_varid, 'description', 'time array');
netcdf.putAtt(ncid, time_varid, 'units', 's');
netcdf.putAtt(ncid, time_varid, 'array dimensions', size(t_arr));

epv_varid = netcdf.defVar(ncid, 'epv_t_hslice', 'double', [t_len y_len x_len]);
netcdf.putAtt(ncid, epv_varid, 'description', 'epv hslice at const vlevel (Nt)');
netcdf.putAtt(ncid, epv_varid, 'units', 's^{-1}');
netcdf.putAtt(ncid, epv_varid, 'array dimensions', size(epv_t_hslice));
%close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, time_varid, t_arr); % start from 0
netcdf.putVar(ncid, epv_varid, epv_t_hslice);
% close netcdf file
netcdf.close(ncid);
disp('Done creating epv_t_hslice at const vlevel!')
%--------------------------------------------------------
%%

   
