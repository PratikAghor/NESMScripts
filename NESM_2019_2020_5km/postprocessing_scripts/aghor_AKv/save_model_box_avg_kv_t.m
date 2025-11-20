%-------------------------------------------------------------------
% save horizontally averaged Kv = <Kv> versus t at a constant vlevel
% Kv [m^2/s] is saved as Akt in the model output
% Author: Pratik Aghor
%-------------------------------------------------------------------
clear all; close all; clc;
start_paths;
%-------------------------------------------------------------------
indxRange = 0:3877; % entire year
nt0 = indxRange(1);
[~, Nt] = size(indxRange);

vlevel = -4000;

domain_avg_kv_t = zeros(Nt, 1);
box_avg_kv_t = zeros(Nt, 1);
t_arr = zeros(Nt, 1);
% nanprofile_counter = 0; % count of nan profiles (due to topography)

box="AtlantisII_5km" % average around AtlantisII
% box="domain_5km" % average around all three mountains, bigger box

[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)

for nt = 1:Nt
        sprintf(strcat('loading', [dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], 'file'), indxRange(1, nt))
        fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
        hisfile = fname;

        t_arr(nt) = ncread(hisfile, 'time');
	
	% AKt = pagetranspose(ncread(hisfile, 'AKt'));
        % AKt = shiftdim(AKt, 2);
        
	% AKv = ncread(hisfile, 'AKv');
	AKv = ncread(hisfile, 'AKt');

	if(nt == 0)
		size(AKv)
	end

	AKv_vlevel = vinterp(AKv, zw, vlevel); % AKv at vlevel
	box_AKv_vlevel = AKv_vlevel(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);
	
	domain_avg_kv_t(nt) = mean(AKv_vlevel, 'all', 'omitnan');
	box_avg_kv_t(nt) = mean(box_AKv_vlevel, 'all', 'omitnan');

	sprintf('domain_avg_kv_t = %.7e', domain_avg_kv_t(nt) )
	sprintf('box_avg_kv_t = %.7e', box_avg_kv_t(nt) )
	
end
%%------------------------------------------
% save timeseries of box avg into a netcdf file
if(strcmp(box, 'AtlantisII_5km')) % covers only Atlantis II
	filename = strcat('nesm_2019_2020_5km_model_box_avg_kv_t', '_vlevel_', string(vlevel), '_nt_', string(indxRange(1)), '_', ...
		string(indxRange(Nt)), '.nc')
elseif(strcmp(box, 'domain_5km')) % covers all three seamounts
	filename = strcat('nesm_2019_2020_5km_model_box_avg_kv_t', '_vlevel_', string(vlevel), '_nt_', string(indxRange(1)), '_', ...
                string(indxRange(Nt)), '.nc')
end
ncid = netcdf.create(filename,'CLOBBER');
%% freq_len = netcdf.defDim(ncid, 'freq_len', length(Fv1));
t_len = netcdf.defDim(ncid, 'Nt', Nt);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
% z_len = netcdf.defDim(ncid, 'Nz', Nz);
netcdf.close(ncid);
%%
% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

t_varid = netcdf.defVar(ncid, 't_arr', 'double', [t_len]);
netcdf.putAtt(ncid, t_varid, 'description', 'time array of time measured from 01/01/2018');
netcdf.putAtt(ncid, t_varid, 'units', 's');
netcdf.putAtt(ncid, t_varid, 'array dimensions', size(t_arr));

domain_avg_k_t_varid = netcdf.defVar(ncid, 'domain_avg_kv_t', 'double', [t_len]);
netcdf.putAtt(ncid, domain_avg_k_t_varid, 'description', 'timeseries of domain averaged Kv near bottom diffusivity');
netcdf.putAtt(ncid, domain_avg_k_t_varid, 'units', 'm^2s^-1');
netcdf.putAtt(ncid, domain_avg_k_t_varid, 'array dimensions', size(domain_avg_kv_t));

k_t_varid = netcdf.defVar(ncid, 'box_avg_kv_t', 'double', [t_len]);
netcdf.putAtt(ncid, k_t_varid, 'description', 'timeseries of box averaged Kv near bottom diffusivity');
netcdf.putAtt(ncid, k_t_varid, 'units', 'm^2s^-1');
netcdf.putAtt(ncid, k_t_varid, 'array dimensions', size(box_avg_kv_t));

% close define mode
netcdf.endDef(ncid);
%%
% put values
netcdf.putVar(ncid, t_varid, t_arr);
netcdf.putVar(ncid, domain_avg_k_t_varid, domain_avg_kv_t);
netcdf.putVar(ncid, k_t_varid, box_avg_kv_t);
% close netcdf file
netcdf.close(ncid);
disp('Done creating near bottom diffusivity!')
%%------------------------------------------
