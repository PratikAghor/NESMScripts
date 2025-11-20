%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create and fill CROCO initial file from hycom outputs.
% by sdxmonkey
% 2020-02-18
%
% modified from make_OGCM of CROCOTOOLS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('/storage/scratch1/9/paghor3/croco/croco_tools')
addpath('/storage/home/hcoda1/9/paghor3/croco/croco_tools')
start;
clear
global nctbx_options;
nctbx_options.theAutoNaN = 1;
nctbx_options.theAutoscale = 1;
%% parameters --------------------------------------------------
CROCO_title  = 'NESM_1km_s100';
CROCO_config = 'NESM_1km_s100';
dir_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/NESM_2019_2020/'
blk_bry_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/blk_bry_data/'

ininame = [dir_path, 'forcing/NESM_ini_20190101_s100.nc']
offset = 0; % offset for the time -1 for p1 case, 1 for m1 case
file_OGCM = [dir_path, 'forcing/hycom_OGCM_20190101_20200501.nc']
Yorig = 2018;
grdname = [dir_path, 'forcing/NESM_grd.nc']
tin = 9;
tout = 1;

% grid information
vtransform=2;
Vstretching=4;
theta_s=3;
theta_b=5;
hc=30; %Tcline=400
N=100;
%
% Objective analysis decorrelation scale [m]
% (if Roa=0: nearest extrapolation method; crude but much cheaper)
%
%Roa=300e3;
Roa=0;
%
interp_method = 'spline';         % Interpolation method: 'linear' or 'cubic'
%end of parapemters

%------------------------------------------------------------------

%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)
[M,L]=size(lon);
%
%
% CROCO grid angle
%
cosa=cos(angle);
sina=sin(angle);

%
% Get the OGCM grid
%
nc=netcdf(file_OGCM);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=nc{'depth'}(:); 
NZ=length(Z);
% NZ=NZ-rmdepth;
Z=Z(1:NZ);
time = nc{'time'}(tin)+offset;
close(nc)

%
% Initial file
%
create_inifile(ininame,grdname,CROCO_title,...
                 theta_s,theta_b,hc,N,...
                 time,'clobber', vtransform);
nc_ini=netcdf(ininame,'write');
% interpolation (modified from interp_OGCM.m) -----------------------
conserv=1; % same barotropic velocities as the OGCM
%
disp(['  Horizontal interpolation: ',...
      file_OGCM])
%
%
% Open the OGCM file
%
nc=netcdf(file_OGCM);
%
% Interpole data on the OGCM Z grid and CROCO horizontal grid
%
% Get zeta because it is needed to compute vertical levels of CROCO grid
zeta=ext_data_OGCM(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method);
%
%
% Read and extrapole the 2D variables
%
  u2d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,lon,lat,1,Roa,interp_method);
  v2d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,lon,lat,1,Roa,interp_method);
  ubar=rho2u_2d(u2d.*cosa+v2d.*sina);
  vbar=rho2v_2d(v2d.*cosa-u2d.*sina);
%
% Read and extrapole the 3D variables
%
  NZ=length(Z);
  dz=gradient(Z);
  temp=zeros(NZ,M,L);
  salt=zeros(NZ,M,L);
  u=zeros(NZ,M,L-1);
  v=zeros(NZ,M-1,L);
  for k=1:NZ
    if rem(k,10)==0
      disp(['  Level ',num2str(k),' of ',num2str(NZ)])
    end
    u2d=ext_data_OGCM(nc,lonU,latU,'u',tin,lon,lat,k,Roa,interp_method);
    v2d=ext_data_OGCM(nc,lonV,latV,'v',tin,lon,lat,k,Roa,interp_method);
    u(k,:,:)=rho2u_2d(u2d.*cosa+v2d.*sina);
    v(k,:,:)=rho2v_2d(v2d.*cosa-u2d.*sina);
    temp(k,:,:)=ext_data_OGCM(nc,lonT,latT,'temp',tin,lon,lat,k,Roa,interp_method);
    salt(k,:,:)=ext_data_OGCM(nc,lonT,latT,'salt',tin,lon,lat,k,Roa,interp_method);
  end
%
% Close the OGCM file
%
close(nc)
%
%
% Get the CROCO vertical grid
%
disp('  Vertical interpolations')
  theta_s=nc_ini{'theta_s'}(:);
  theta_b=nc_ini{'theta_b'}(:);
  hc=nc_ini{'hc'}(:);
  N=length(nc_ini('s_rho'));
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
Z=[100;Z;-100000];
%
% CROCO vertical grid
%

zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
zu=rho2u_3d(zr);
zv=rho2v_3d(zr);
zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=rho2u_3d(dzr);
dzv=rho2v_3d(dzr);

%
% Add a level on top and bottom with no-gradient
%
  u=cat(1,u(1,:,:),u);
  u=cat(1,u,u(end,:,:));
  v=cat(1,v(1,:,:),v);
  v=cat(1,v,v(end,:,:));
  temp=cat(1,temp(1,:,:),temp);
  temp=cat(1,temp,temp(end,:,:));
  salt=cat(1,salt,salt(end,:,:));
  salt=cat(1,salt(1,:,:),salt);

%
% Perform the vertical interpolations
% aghoredit: used ztosigma_1d instead of ztosigma
%
  u=ztosigma(flipdim(u,1),zu,flipud(Z));
  v=ztosigma(flipdim(v,1),zv,flipud(Z));
  temp=ztosigma(flipdim(temp,1),zr,flipud(Z));
  salt=ztosigma(flipdim(salt,1),zr,flipud(Z));
%
% Correct the horizontal transport
% i.e. remove the interpolated tranport and add
%      the OGCM transport
%
  if conserv==1
    u=u-tridim(squeeze(sum(u.*dzu)./sum(dzu)),N);
    v=v-tridim(squeeze(sum(v.*dzv)./sum(dzv)),N);
    u=u+tridim(ubar,N);
    v=v+tridim(vbar,N);
  end
%
% Barotropic velocities
%
  ubar=squeeze(sum(u.*dzu)./sum(dzu));
  vbar=squeeze(sum(v.*dzv)./sum(dzv));
%
%  fill the files
%
  nc_ini{'zeta'}(tout,:,:)=zeta;
  nc_ini{'SSH'}(tout,:,:)=zeta;
  nc_ini{'temp'}(tout,:,:,:)=temp;
  nc_ini{'salt'}(tout,:,:,:)=salt;
  nc_ini{'u'}(tout,:,:,:)=u;
  nc_ini{'v'}(tout,:,:,:)=v;
  nc_ini{'ubar'}(tout,:,:,:)=ubar;
  nc_ini{'vbar'}(tout,:,:,:)=vbar;
%--------------------------------------------------------------------
close(nc_ini)

disp('done creating inifile');
