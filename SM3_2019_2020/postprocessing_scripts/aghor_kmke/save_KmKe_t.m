% 
% Compute the kinetic energy transfer: save KmKe3D time series at a given depth
% - Pratik Aghor, modified from get_KmKe.m
% HRS= -[ <up up> dubar/dx + <up vp>  dubar/dy  ....
%          <vp up> dvbar/dx + <vp vp>  dvbar/dy ]
%     =  -<up(up.grad(ubar))+vp(up.grad(vbar))>
%
% Advection operators are used. Much easier in sigma coordinates.
%
% This is the method used in
% Djakouré, S., P. Penven, B. Bourlès, J. Veitch and V. Koné, 
% Coastally trapped eddies in the north of the Gulf of Guinea, 2014, J. Geophys. Res,
% DOI: 10.1002/2014JC010243
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Assume annual mean has been calculated and saved
clear all
close all
start_paths
%-------------------------------------------------------
% box limits for horizontal averages
box="AtlantisII_1km" % average around AtlantisII
% box="NESM_1km"
% box = "domain_1km"
[lat_idx_min, lat_idx_max, lon_idx_min, lon_idx_max] = box_idx_lims(box)
box_Ny = lat_idx_max - lat_idx_min
box_Nx = lon_idx_max - lon_idx_min
%-------------------------------------------------------

%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;

vlevel = -3000;

% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % entire year
nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
ub=zeros(N, Mu, Lu);
vb=zeros(N, Mv, Lv);

annual_mean_indxRange = 952:3877;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_uvw_file = strcat('sm3_2019_2020_uvw_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

ub = ncread(annual_mean_uvw_file, 'ub');
vb = ncread(annual_mean_uvw_file, 'vb');
wb = ncread(annual_mean_uvw_file, 'wb');
zetab = ncread(annual_mean_uvw_file, 'zetab');

[Wvlcb,Wb]=get_wvelocity(zetab,ub,vb,depth,pm,pn,theta_s,theta_b,hc,N,Vtransform);

%%
%
% gradients of ubar and vbar
%

dudx=zeros(N,M,L);
dudx(:,:,2:end-1)=tridim(pm(:,2:end-1),N).*(ub(:,:,2:end)-ub(:,:,1:end-1));
dudy=zeros(N,M,L);
cff=ub(:,2:end,:)-ub(:,1:end-1,:);
dudy(:,2:end-1,2:end-1)=tridim(pn(2:end-1,2:end-1),N).*0.25.*...
             (cff(:,2:end,2:end)  +cff(:,1:end-1,2:end)+...
              cff(:,2:end,1:end-1)+cff(:,1:end-1,1:end-1));
dvdx=zeros(N,M,L);
cff=vb(:,:,2:end)-vb(:,:,1:end-1);
dvdx(:,2:end-1,2:end-1)=tridim(pm(2:end-1,2:end-1),N).*0.25.*...
             (cff(:,2:end,2:end)  +cff(:,1:end-1,2:end)+...
              cff(:,2:end,1:end-1)+cff(:,1:end-1,1:end-1));
dvdy=zeros(N,M,L);
dvdy(:,2:end-1,:)=tridim(pn(2:end-1,:),N).*(vb(:,2:end,:)-vb(:,1:end-1,:));
%
dn_u=tridim(2./(pn(:,1:end-1)+pn(:,2:end)),N);
dm_v=tridim(2./(pm(1:end-1,:)+pm(2:end,:)),N);
omn_w=tridim(1./(pm.*pn),N+1);

%
% Time loop
%
% KmKe_t = zeros(Nt, N, M, L);
KmKe_t_box_avg = zeros(Nt, 1);
HRS_t_box_avg = zeros(Nt, 1);
VRS_t_box_avg = zeros(Nt, 1);

for nt=1:Nt
    HRS=zeros(N,M,L);
    VRS=zeros(N,M,L);

    fname = sprintf([dirHR, 'SM3_2019_2020_avg.%05d.nc'], indxRange(1, nt));
    sprintf(strcat('loading', fname))
    hisfile = fname;
    % nc=netcdf(fname);
    % u=squeeze(nc{'u'}(nt,:,:,:));
    % v=squeeze(nc{'v'}(nt,:,:,:));
    % calculate HRS and VRS
    u = pagetranspose(ncread(hisfile, 'u'));
    v = pagetranspose(ncread(hisfile, 'v'));
    u = shiftdim(u, 2);
    v = shiftdim(v, 2);
    zeta = pagetranspose(ncread(hisfile, 'zeta'));
    up = u - ub; % u' = u - \bar{u}
    vp = v - vb; % v' = v - \bar{v}
  
    %
    % Compute Hz
    % 
    zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);
    Hz=zw(2:N+1,:,:)-zw(1:N,:,:); 
    mnoHz=tridim(pm.*pn,N)./Hz;
    %
    % Compute horizontal advection terms
    %
    trc_u = u2rho_3d(ub);
    % shift dimensions to pass into trcflux fun
    % trc_u = shiftdim(trc_u, 1);
    % trc_u = pagetranspose(trc_u);
    
    % size(trc_u)
    upXgradub=mnoHz.*croco_horiz_advection_aghor(masku, maskv, maskr, ...
	Hz,dn_u.*up, dm_v.*vp, trc_u);

    trc_v = v2rho_3d(vb);
    % shift dims for compatibility with trcflux
    % trc_v = shiftdim(trc_v, 1);
    % trc_v = pagetranspose(trc_v);
    upXgradvb=mnoHz.*croco_horiz_advection_aghor(masku, maskv, maskr, ...
	Hz,dn_u.*up,dm_v.*vp,trc_v);

    % already calculated at rho pts, no need to convert
    up=u2rho_3d(up);
    vp=v2rho_3d(vp);

    HRS= - up.*upXgradub - vp.*upXgradvb;
    %
    % VRS
    [Wvlc,Wrk]=get_wvelocity(zeta,u,v,depth,pm,pn,theta_s,theta_b,hc,N,Vtransform);


    size(Wrk)
    Wp=Wrk-Wb;
    %
    % Compute Hz, already computed above
    % 
    %
    % Compute vertical advection terms
    %
    wpXgradub=mnoHz.*croco_vert_advection_aghor(maskr, omn_w.*Wp, u2rho_3d(ub));
    wpXgradvb=mnoHz.*croco_vert_advection_aghor(maskr, omn_w.*Wp, v2rho_3d(vb));

    % up, vp already computed above
    % up=u2rho_3d(up);
    % vp=v2rho_3d(vp);

    VRS = - up.*wpXgradub - vp.*wpXgradvb;
    %
    % KmKe_t(nt, :, :, :) = HRS + VRS;
    KmKe_temp = HRS + VRS;

    % updated from save_KmKe_2d.m, take an hslice at a given vlevel
    zeta = pagetranspose(ncread(hisfile, 'zeta'));

    zr=zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',Vtransform);
    zw=zlevs(depth,zeta,theta_s,theta_b,hc,N,'w',Vtransform);

    KmKe2D = vinterp(KmKe_temp,zr,vlevel);
    HRS2D = vinterp(HRS, zr, vlevel);
    VRS2D = vinterp(VRS, zr, vlevel);

    % take a horizontal avg around the middle seamount
    KmKe_t_box = KmKe2D(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);
    HRS_t_box = HRS2D(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);
    VRS_t_box = VRS2D(lat_idx_min:lat_idx_max, lon_idx_min:lon_idx_max);

    KmKe_t_box_avg(nt) = mean(KmKe_t_box, 'all', 'omitnan');
    HRS_t_box_avg(nt) = mean(HRS_t_box, 'all', 'omitnan');
    VRS_t_box_avg(nt) = mean(VRS_t_box, 'all', 'omitnan');

end

%
% Save
%
% save data into a netcdf file
filename = strcat('sm3_2019_2020_KmKe_t_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '_vlevel_', string(vlevel), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
t_len = netcdf.defDim(ncid, 'Nt', Nt);
z_len = netcdf.defDim(ncid, 'Nz', N);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
x_len = netcdf.defDim(ncid, 'Nx', Nx);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

hrs_varid = netcdf.defVar(ncid, 'HRS_t', 'double', [t_len]);
netcdf.putAtt(ncid, hrs_varid, 'description', 'box averaged HRS timeseries (Nt) at a const depth');
netcdf.putAtt(ncid, hrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, hrs_varid, 'array dimensions', size(HRS_t_box_avg));

vrs_varid = netcdf.defVar(ncid, 'VRS_t', 'double', [t_len]);
netcdf.putAtt(ncid, vrs_varid, 'description', 'box averaged VRS timeseries (Nt) at a const depth');
netcdf.putAtt(ncid, vrs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, vrs_varid, 'array dimensions', size(VRS_t_box_avg));


rs_varid = netcdf.defVar(ncid, 'KmKe_t', 'double', [t_len]);
netcdf.putAtt(ncid, rs_varid, 'description', 'box averaged KmKe timeseries (Nt) at a const depth');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, rs_varid, 'array dimensions', size(KmKe_t_box_avg));
% close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, z01_varid, z01); % start from 0
% netcdf.putVar(ncid, z02_varid, z02); % start from 0
netcdf.putVar(ncid, hrs_varid, HRS_t_box_avg); % start from 0
netcdf.putVar(ncid, vrs_varid, VRS_t_box_avg); % start from 0
netcdf.putVar(ncid, rs_varid, KmKe_t_box_avg); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating KmKe_t!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


