%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create and fill CROCO bry file from hycom outputs.
% by sdxmonkey
% 2020-02-18
%
% modified from make_OGCM of CROCOTOOLS
%
% better run with nohup:
% nohup /home/PublicSoftware/matlab-R2017B/bin/matlab -nodesktop -nosplash -nodisplay < make_BRY_hycom.m > make_BRY.log 2>&1 &
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('/storage/scratch1/9/paghor3/croco/croco_tools/');
addpath('/storage/home/hcoda1/9/paghor3/croco/croco_tools/');
start
clear
global nctbx_options;
nctbx_options.theAutoNaN = 1;
nctbx_options.theAutoscale = 1;
%% parameters --------------------------------------------------
CROCO_title  = 'NESM_1km_s100';
CROCO_config = 'NESM_1km_s100';

dir_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/NESM_2019_2020/'
blk_bry_path = '/storage/home/hcoda1/9/paghor3/scratch/work/new_runs_2/blk_bry_data/'

bryname = [dir_path, 'forcing/NESM_bry_20190101_20200501_s100.nc'];
file_OGCM = [dir_path, 'forcing/hycom_OGCM_20190101_20200501.nc']; 
Yorig = 2018;
grdname = [dir_path, 'forcing/NESM_grd.nc']
ifcreate = 1; % 0 once OGCM file is written
file_orig_begin = [blk_bry_path, 'merged/'];
file_orig_end = '_20190101_20200501.nc'
obc = [1 1 1 1] % open boundaries (1=open , [S E N W])
obc_names = {'south','east','north','west'};
IIobc = find(obc == 1);
t0 = 0;
Nt = 3884; % nrecords

% grid information
vtransform=2;
%Vstretching=1;
theta_s=3;
theta_b=5;
hc=30; %Tcline=400
N=100;

%
% Objective analysis decorrelation scale [m]
% (if Roa=0: nearest extrapolation method; crude but much cheaper)
%
%Roa=300e3;
%Roa=0;
Roa = 10e3;
%
interp_method = 'spline';         % Interpolation method: 'linear' or 'cubic'
%end of parapemters

%------------------------------------------------------------------

% preprocess original hycom file
if ifcreate
    disp(['Preprocessing ...'])
    create_hycom_bry_all(file_orig_begin,file_orig_end,file_OGCM,Yorig);
end

%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
% CROCO grid angle
cosa=cos(angle);
sina=sin(angle);
close(nc)
[M,L]=size(lon);

nc=netcdf(file_OGCM);
time = nc{'time'}(:);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=nc{'depth'}(:);
NZ=length(Z);
%NZ=NZ-rmdepth;
Z=Z(1:NZ);
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
Z=[100;Z;-100000];
%
% Initialize the file
%
create_bryfile(bryname,grdname,CROCO_title,obc,...
                       theta_s,theta_b,hc,N,...
                       time(t0+1:t0+Nt),0,'clobber',vtransform);
nc_bry=netcdf(bryname,'write');
ntimes = length(time);

%
% Get the CROCO vertical grid
%
theta_s=nc_bry{'theta_s'}(:);
theta_b=nc_bry{'theta_b'}(:);
hc=nc_bry{'hc'}(:);
N=length(nc_bry('s_rho'));

for tndx_OGCM = 1:Nt
  disp(['Time step : ',num2str(tndx_OGCM),' of ',num2str(Nt),' :'])
  tin = tndx_OGCM+t0;
  tout = tndx_OGCM;
  conserv=1; % same barotropic velocities as the OGCM
%
  for ii = 1:length(IIobc)
    disp(['  Working on: ',obc_names(IIobc(ii))])
    %
    % Get the OGCM grid
    %
    disp(['    Horizontal interpolation ...'])
   %
    % Interpole data on the OGCM Z grid and CROCO horizontal grid
    %
    % Get zeta because it is needed to compute vertical levels of CROCO grid
    zeta=ext_data_OGCM(nc,lonT,latT,'ssh',tin,lon,lat,1,Roa,interp_method);
    if IIobc(ii) == 1
      zeta_south=squeeze(zeta(1,:));
      ubar1d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,...
           squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),1,Roa,interp_method);
      vbar1d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,...
           squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),1,Roa,interp_method);
      ubar_south=squeeze(rho2u_2d(ubar1d(1,:).*cosa(1,:)+vbar1d(1,:).*sina(1,:)));
      vbar_south=squeeze(rho2v_2d(vbar1d.*cosa(1:2,:)-ubar1d.*sina(1:2,:)));
    end
    if IIobc(ii) == 2
      zeta_east=squeeze(zeta(:,end));
      ubar1d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,...
           squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),1,Roa,interp_method);
      vbar1d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,...
           squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),1,Roa,interp_method);
      ubar_east=squeeze(rho2u_2d(ubar1d.*cosa(:,end-1:end)+vbar1d.*sina(:,end-1:end)));
      vbar_east=squeeze(rho2v_2d(vbar1d(:,end).*cosa(:,end)-ubar1d(:,end).*sina(:,end)));
    end
    if IIobc(ii) == 3
      zeta_north=squeeze(zeta(end,:));
      ubar1d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,...
           squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),1,Roa,interp_method);
      vbar1d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,...
           squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),1,Roa,interp_method);
      ubar_north=squeeze(rho2u_2d(ubar1d(end,:).*cosa(end,:)+vbar1d(end,:).*sina(end,:)));
      vbar_north=squeeze(rho2v_2d(vbar1d.*cosa(end-1:end,:)-ubar1d.*sina(end-1:end,:)));
    end
    if IIobc(ii) == 4
      zeta_west=squeeze(zeta(:,1));
      ubar1d=ext_data_OGCM(nc,lonU,latU,'ubar',tin,...
           squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),1,Roa,interp_method);
      vbar1d=ext_data_OGCM(nc,lonV,latV,'vbar',tin,...
           squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),1,Roa,interp_method);
      ubar_west=squeeze(rho2u_2d(ubar1d.*cosa(:,1:2)+vbar1d.*sina(:,1:2)));
      vbar_west=squeeze(rho2v_2d(vbar1d(:,1).*cosa(:,1)-ubar1d(:,1).*sina(:,1)));
    end
    %
    % Read and extrapole the 3D variables
    %
    if IIobc(ii) == 1
      temp_south=zeros(NZ,L);
      salt_south=zeros(NZ,L);
      u_south=zeros(NZ,L-1);
      v_south=zeros(NZ,L);
    end
    if IIobc(ii) == 2
      temp_east=zeros(NZ,M);
      salt_east=zeros(NZ,M);
      u_east=zeros(NZ,M);
      v_east=zeros(NZ,M-1);
    end
    if IIobc(ii) == 3
      temp_north=zeros(NZ,L);
      salt_north=zeros(NZ,L);
      u_north=zeros(NZ,L-1);
      v_north=zeros(NZ,L);
    end  
    if IIobc(ii) == 4
      temp_west=zeros(NZ,M);
      salt_west=zeros(NZ,M);
      u_west=zeros(NZ,M);
      v_west=zeros(NZ,M-1);
    end
  %
    for k=1:NZ
      if rem(k,10)==0
        disp(['  Level bry ',num2str(k),' of ',num2str(NZ)])
      end
      if  IIobc(ii) == 1 % Southern boundary
        t1d=squeeze(ext_data_OGCM(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k,Roa,interp_method));
        s1d=squeeze(ext_data_OGCM(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k,Roa,interp_method));
        temp_south(k,:)=squeeze(t1d(1,:));
        salt_south(k,:)=squeeze(s1d(1,:));
        u1d=ext_data_OGCM(nc,lonU,latU,'u',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k,Roa,interp_method);
        v1d=ext_data_OGCM(nc,lonV,latV,'v',tin,...
                  squeeze(lon(1:2,:)),squeeze(lat(1:2,:)),...
                  k,Roa,interp_method);
        u_south(k,:)=squeeze(rho2u_2d(u1d(1,:).*cosa(1,:)+...
                                    v1d(1,:).*sina(1,:)));
        v_south(k,:)=squeeze(rho2v_2d(v1d.*cosa(1:2,:)-...
                                    u1d.*sina(1:2,:)));
      end
      if IIobc(ii) == 2  % Eastern boundary
        t1d=squeeze(ext_data_OGCM(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k,Roa,interp_method));
        s1d=squeeze(ext_data_OGCM(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k,Roa,interp_method));
        temp_east(k,:)=squeeze(t1d(:,end)');
        salt_east(k,:)=squeeze(s1d(:,end)');
        u1d=ext_data_OGCM(nc,lonU,latU,'u',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k,Roa,interp_method);
        v1d=ext_data_OGCM(nc,lonV,latV,'v',tin,...
                  squeeze(lon(:,end-1:end)),squeeze(lat(:,end-1:end)),...
                  k,Roa,interp_method);
        cff=squeeze(rho2u_2d(u1d.*cosa(:,end-1:end)+...
                           v1d.*sina(:,end-1:end)));
        u_east(k,:)=cff';
        cff=squeeze(rho2v_2d(v1d(:,end).*cosa(:,end)-...
                           u1d(:,end).*sina(:,end)));
        v_east(k,:)=cff';
      end
      if IIobc(ii) == 3  % Northern boundary
        t1d=squeeze(ext_data_OGCM(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k,Roa,interp_method));
        s1d=squeeze(ext_data_OGCM(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k,Roa,interp_method));
        temp_north(k,:)=squeeze(t1d(end,:));
        salt_north(k,:)=squeeze(s1d(end,:));
        u1d=ext_data_OGCM(nc,lonU,latU,'u',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k,Roa,interp_method);
        v1d=ext_data_OGCM(nc,lonV,latV,'v',tin,...
                  squeeze(lon(end-1:end,:)),squeeze(lat(end-1:end,:)),...
                  k,Roa,interp_method);
        u_north(k,:)=squeeze(rho2u_2d(u1d(end,:).*cosa(end,:)+...
                                    v1d(end,:).*sina(end,:)));
        v_north(k,:)=squeeze(rho2v_2d(v1d.*cosa(end-1:end,:)-...
                                    u1d.*sina(end-1:end,:)));
      end
      if IIobc(ii) == 4  % Western boundary
        t1d=squeeze(ext_data_OGCM(nc,lonT,latT,'temp',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k,Roa,interp_method));
        s1d=squeeze(ext_data_OGCM(nc,lonT,latT,'salt',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k,Roa,interp_method));
        temp_west(k,:)=squeeze(t1d(:,1)');
        salt_west(k,:)=squeeze(s1d(:,1)');
        u1d=ext_data_OGCM(nc,lonU,latU,'u',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k,Roa,interp_method);
        v1d=ext_data_OGCM(nc,lonV,latV,'v',tin,...
                  squeeze(lon(:,1:2)),squeeze(lat(:,1:2)),...
                  k,Roa,interp_method);
        cff=squeeze(rho2u_2d(u1d.*cosa(:,1:2)+...
                           v1d.*sina(:,1:2)));
        u_west(k,:)=cff';
        cff=squeeze(rho2v_2d(v1d(:,1).*cosa(:,1)-...
                           u1d(:,1).*sina(:,1)));
        v_west(k,:)=cff';
      end
    end
    %
    % CROCO vertical grid
    %
    % aghoredit: only needed when h is negative in grd
    % zlevs has Dcrit - h which makes zr positive to give
    % erros in the vertical interpilations
    % therefore, just to get zr have negative vals:
    % h = -h;
    zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);
    zu=rho2u_3d(zr);
    zv=rho2v_3d(zr);
    zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',vtransform);
    dzr=zw(2:end,:,:)-zw(1:end-1,:,:);
    dzu=rho2u_3d(dzr);
    dzv=rho2v_3d(dzr);
    % aghoredit, reset h to original vals: only needed when h is negative in grd
    % h = -h;
    %
    % Vertical interpolation in case of bry files
    %
    if IIobc(ii) == 1
      [u_south,v_south,ubar_south,vbar_south,...
       temp_south,salt_south]=vinterp_OGCM_bry(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
                                             dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
                                             u_south,v_south,ubar_south,vbar_south,...
                                             temp_south,salt_south,...
                                             N,Z,conserv);
    end
    if IIobc(ii) == 2
      [u_east,v_east,ubar_east,vbar_east,...
       temp_east,salt_east]=vinterp_OGCM_bry(zr(:,:,end),zu(:,:,end),zv(:,:,end),...
                                           dzr(:,:,end),dzu(:,:,end),dzv(:,:,end),...
                                           u_east,v_east,ubar_east,vbar_east,...
                                           temp_east,salt_east,...
                                           N,Z,conserv);
    end
    if IIobc(ii) == 3
      [u_north,v_north,ubar_north,vbar_north,...
       temp_north,salt_north]=vinterp_OGCM_bry(zr(:,end,:),zu(:,end,:),zv(:,end,:),...
                                             dzr(:,end,:),dzu(:,end,:),dzv(:,end,:),...
                                             u_north,v_north,ubar_north,vbar_north,...
                                             temp_north,salt_north,...
                                             N,Z,conserv);
    end
    if IIobc(ii) == 4
      [u_west,v_west,ubar_west,vbar_west,...
       temp_west,salt_west]=vinterp_OGCM_bry(zr(:,:,1),zu(:,:,1),zv(:,:,1),...
                                           dzr(:,:,1),dzu(:,:,1),dzv(:,:,1),...
                                           u_west,v_west,ubar_west,vbar_west,...
                                           temp_west,salt_west,...
                                           N,Z,conserv);
    end
    %
    %  fill the files
    %
    if IIobc(ii) == 1
      nc_bry{'zeta_south'}(tout,:)=zeta_south;
      nc_bry{'temp_south'}(tout,:,:)=temp_south;
      nc_bry{'salt_south'}(tout,:,:)=salt_south;
      nc_bry{'u_south'}(tout,:,:)=u_south;
      nc_bry{'v_south'}(tout,:,:)=v_south;
      nc_bry{'ubar_south'}(tout,:,:)=ubar_south;
      nc_bry{'vbar_south'}(tout,:,:)=vbar_south;
    end 
    if IIobc(ii) == 2
      nc_bry{'zeta_east'}(tout,:)=zeta_east;
      nc_bry{'temp_east'}(tout,:,:)=temp_east;
      nc_bry{'salt_east'}(tout,:,:)=salt_east;
      nc_bry{'u_east'}(tout,:,:)=u_east;
      nc_bry{'v_east'}(tout,:,:)=v_east;
      nc_bry{'ubar_east'}(tout,:,:)=ubar_east;
      nc_bry{'vbar_east'}(tout,:,:)=vbar_east;
    end 
    if IIobc(ii) == 3
      nc_bry{'zeta_north'}(tout,:)=zeta_north;
      nc_bry{'temp_north'}(tout,:,:)=temp_north;
      nc_bry{'salt_north'}(tout,:,:)=salt_north;
      nc_bry{'u_north'}(tout,:,:)=u_north;
      nc_bry{'v_north'}(tout,:,:)=v_north;
      nc_bry{'ubar_north'}(tout,:,:)=ubar_north;
      nc_bry{'vbar_north'}(tout,:,:)=vbar_north;
    end 
    if IIobc(ii) == 4
      nc_bry{'zeta_west'}(tout,:)=zeta_west;
      nc_bry{'temp_west'}(tout,:,:)=temp_west;
      nc_bry{'salt_west'}(tout,:,:)=salt_west;
      nc_bry{'u_west'}(tout,:,:)=u_west;
      nc_bry{'v_west'}(tout,:,:)=v_west;
      nc_bry{'ubar_west'}(tout,:,:)=ubar_west;
      nc_bry{'vbar_west'}(tout,:,:)=vbar_west;
    end
  end
end
close(nc)
close(nc_bry)
    %------------------------------------------------------------------------------
