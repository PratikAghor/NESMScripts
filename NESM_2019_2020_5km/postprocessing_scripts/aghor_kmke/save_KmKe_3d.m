% 
% Compute the kinetic energy transfer: save KmKe timeseries (not time average) at a given depth
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
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
% indxRange = 469:3382; % What time indices do you need?
% indxRange = 2120:2200;
indxRange = 2984:3040;
nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
ub=zeros(N, Mu, Lu);
vb=zeros(N, Mv, Lv);

annual_mean_indxRange = 469:3382;
[~, mean_Nt] = size(annual_mean_indxRange);
annual_mean_uvw_file = strcat('nesm_2019_2020_uvw_annual_mean_3d_nt_',string(annual_mean_indxRange(1)), '_', string(annual_mean_indxRange(mean_Nt)),'.nc');

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
KmKe_t = zeros(Nt, N, M, L);
for nt=1:Nt
    HRS=zeros(N,M,L);
    VRS=zeros(N,M,L);

    fname = sprintf([dirHR, 'NESM_2019_2020_avg.%05d.nc'], indxRange(1, nt));
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
    KmKe_t(nt, :, :, :) = HRS + VRS;
    nt=nt+1;
    %
end

%
% Save
%
% save data into a netcdf file
filename = strcat('nesm_2019_2020_KmKe_3d_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
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

rs_varid = netcdf.defVar(ncid, 'KmKe', 'double', [t_len z_len y_len x_len]);
netcdf.putAtt(ncid, rs_varid, 'description', '3d KmKe timeseries (Nt, Nz, Ny, Nx)');
netcdf.putAtt(ncid, rs_varid, 'units', 'm^2s^-3');
netcdf.putAtt(ncid, rs_varid, 'array dimensions', size(KmKe_t));
% close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, z01_varid, z01); % start from 0
% netcdf.putVar(ncid, z02_varid, z02); % start from 0
netcdf.putVar(ncid, rs_varid, KmKe_t); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating KmKe!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


