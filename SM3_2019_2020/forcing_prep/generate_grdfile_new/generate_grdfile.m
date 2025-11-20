%% generate grid file for GOM runs
clear all
clc
% on prometheus
start_paths

%% PART1 --- configurations

%  ROMS title names and directories
ROMS_title  = 'SM3_grd';
NESM_grdfile = [grid_path, 'NESM_grd.nc'];

% Grid dimensions: 
latmin = 37.5;   % Minimum latitudeF  [degree north]
latmax = 40.0;   % Maximum latitude  [degree north]
lonmin = -64.99;   % Minimum longitude [degree east]
lonmax = -61.52;   % Maximum longitude [degree east]
%-------------------------
% % aghoredit
% % Atlantis II lat-lon lims
% atl_lonmin = -63.5;   % Minimum longitude of Atlantis II  [degree east]
% atl_lonmax = -62.84;   % Maximum longitude of Atlantis II [degree east]
% atl_latmin = 38;   % Minimum latitude of Atlantis II [degree north]
% atl_latmax = 39;   % Maximum latitude of Atlantis II [degree north]

%--------------------------
% Grid resolution [degree]
dl = 0.01160; % for 1 km resolution
% dl = 0.0577; % for 5 km resolution
% Number of vertical levels
N = 100;

%  Vertical grid parameters 
theta_s    =  3.;
theta_b    =  5;
hc         = 30.;
vtransform =  2; % s-coordinate type (1: old- ; 2: new- coordinates)

% Minimum depth at the shore [m] (depends on the resolution,
% rule of thumb: dl=1, hmin=300, dl=1/4, hmin=150, ...)
% This affect the filtering since it works on grad(h)/h.
hmin = 30;

% Maximum depth at the shore [m] (to prevent the generation of too big 
% walls along the coast)
hmax_coast = 50;

% Maximum depth [m] (cut the topography to prevent extrapolations below 
% WOA data)
hmax = 5500;

% Slope parameter (r=grad(h)/h) maximum value for topography smoothing
rtarget = 0.25;

% Number of pass of a selective filter to reduce the isolated
% seamounts on the deep ocean.
% n_filter_deep_topo=87;
n_filter_deep_topo=0;

% Number of pass of a single hanning filter at the end of the smooting 
% procedure to ensure that there is no 2DX noise in the topography.
n_filter_final=4;

%  GSHSS user defined coastline (see m_map) 
%  XXX_f.mat    Full resolution data
%  XXX_h.mat    High resolution data
%  XXX_i.mat    Intermediate resolution data
%  XXX_l.mat    Low resolution data
%  XXX_c.mat    Crude resolution data

% coastfileplot = './coastline_i.mat';
% coastfilemask = './coastline_i_mask.mat';

% Objective analysis decorrelation scale [m]
% (if Roa=0: nearest extrapolation method; crude but much cheaper)
Roa=300e3;

interp_method = 'cubic';       % Interpolation method: 'linear' or 'cubic'
makeplot = 1;     % 1: create a few graphics after each preprocessing step

% Set directory
Run_dir = './';
Data_dir = '../../data/'; 
topofile = [Data_dir,'Topo/ETOPO2v2g_f4.nc']; 
% Set grid filename
grdname = [grid_path, ROMS_title,'.nc'];

%% PART2 --- end user configurations
isoctave=exist('octave_config_info');
disp(' ')
disp([' Making the grid: ',grdname])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])
%---------------------------------
% % Get the longitude
% lonr=(lonmin:dl:lonmax);
% 
% % Get the latitude for an isotropic grid
% i=1;
% latr(i)=latmin;
% while latr(i)<=latmax
%   i=i+1;
%   latr(i)=latr(i-1)+dl*cos(latr(i-1)*pi/180);
% end
%--------------------------------
% aghoredit
lonr = ncread(NESM_grdfile, 'lon_rho');
lonr = squeeze(lonr(:, 1));
latr = ncread(NESM_grdfile, 'lat_rho');
latr = squeeze(latr(1, :));
%--------------------------------
[Lonr,Latr]=meshgrid(lonr,latr);
[Lonu,Lonv,Lonp]=rho2uvp(Lonr); 
[Latu,Latv,Latp]=rho2uvp(Latr);

% Create the grid file
disp(' ')
disp(' Create the grid file...')
[M,L]=size(Latp);
disp([' LLm = ',num2str(L-1)])
disp([' MMm = ',num2str(M-1)])
create_grid(L,M,grdname,ROMS_title)

% Fill the grid file
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'lat_u'}(:)=Latu;
nc{'lon_u'}(:)=Lonu;
nc{'lat_v'}(:)=Latv;
nc{'lon_v'}(:)=Lonv;
nc{'lat_rho'}(:)=Latr;
nc{'lon_rho'}(:)=Lonr;
nc{'lat_psi'}(:)=Latp;
nc{'lon_psi'}(:)=Lonp;

% aghoredit
nc{'Vtransform'}(:)=vtransform;
nc{'theta_s'}(:)=theta_s;
nc{'theta_b'}(:)=theta_b;

close(nc)

%  Compute the metrics
disp(' ')
disp(' Compute the metrics...')
[pm,pn,dndx,dmde]=get_metrics(grdname);
xr=0.*pm;
yr=xr;
for i=1:L
  xr(:,i+1)=xr(:,i)+2./(pm(:,i+1)+pm(:,i));
end
for j=1:M
  yr(j+1,:)=yr(j,:)+2./(pn(j+1,:)+pn(j,:));
end
[xu,xv,xp]=rho2uvp(xr);
[yu,yv,yp]=rho2uvp(yr);
dx=1./pm;
dy=1./pn;
dxmax=max(max(dx/1000));
dxmin=min(min(dx/1000));
dymax=max(max(dy/1000));
dymin=min(min(dy/1000));
disp(' ')
disp([' Min dx=',num2str(dxmin),' km - Max dx=',num2str(dxmax),' km'])
disp([' Min dy=',num2str(dymin),' km - Max dy=',num2str(dymax),' km'])

%  Angle between XI-axis and the direction
%  to the EAST at RHO-points [radians].
angle=get_angle(Latu,Lonu);

%  Coriolis parameter
f=4*pi*sin(pi*Latr/180)*366.25/(24*3600*365.25);

% Fill the grid file
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'xl'}(:)=xr(1,end);
nc{'el'}(:)=yr(end,1);
nc{'pm'}(:)=pm;
nc{'pn'}(:)=pn;
nc{'dndx'}(:)=dndx;
nc{'dmde'}(:)=dmde;
nc{'x_u'}(:)=xu;
nc{'y_u'}(:)=yu;
nc{'x_v'}(:)=xv;
nc{'y_v'}(:)=yv;
nc{'x_rho'}(:)=xr;
nc{'y_rho'}(:)=yr;
nc{'x_psi'}(:)=xp;
nc{'y_psi'}(:)=yp;
nc{'angle'}(:)=angle;
nc{'f'}(:)=f;
nc{'spherical'}(:)='T';
close(nc);

%  Add topography from topofile
disp(' ')
disp(' Add topography...')
h=add_topo(grdname,topofile);
%-------------------------------------------
% aghoredit

% 
% h_nesm = ncread(NESM_grdfile, 'h');
% 
h_new=hmax*ones(size(h));
% 
% h_new = generate_gaussian_seamount(latr, lonr, hmax, h);
h_new = generate_gaussian_seamount_chain(latr, lonr, hmax, h);

disp(size(h_new));
disp(size(h));
%-------------------------------------------
% Compute the mask
maskr=h>0;
maskr=process_mask(maskr);
[masku,maskv,maskp]=uvp_mask(maskr);
%  Write it down
nc=netcdf(grdname,'write');
nc{'h'}(:)=h_new;
nc{'mask_u'}(:)=masku;
nc{'mask_v'}(:)=maskv;
nc{'mask_psi'}(:)=maskp;
nc{'mask_rho'}(:)=maskr;
close(nc);

if (isoctave == 0)
% Create the coastline
% if ~isempty(coastfileplot)
%   make_coast(grdname,coastfileplot);
% end

r=input('Do you want to use editmask ? y,[n]','s');
if strcmp(r,'y')
  disp(' Editmask:')
  disp(' Edit manually the land mask.')
  disp(' Press enter when finished.')
  disp(' ')
  if ~isempty(coastfileplot)
    editmask(grdname,coastfilemask)
  else
    editmask(grdname)
  end
  r=input(' Finished with edit mask ? [press enter when finished]','s');
end

close all
end % isoctave

%  Smooth the topography
nc=netcdf(grdname,'write');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);

%  Write it down
disp(' ')
disp(' Write it down...')
nc{'h'}(:)=h;
close(nc);

%% PART3 --- plot figures

% if (isoctave == 0)
% if makeplot==1
%   disp(' ')
%   disp(' Do a plot...')
%   themask=ones(size(maskr));
%   themask(maskr==0)=NaN; 
%   domaxis=[min(min(Lonr)) max(max(Lonr)) min(min(Latr)) max(max(Latr))];
%   colaxis=[min(min(h)) max(max(h))];
%   fixcolorbar([0.25 0.05 0.5 0.03],colaxis,...
%               'Topography',10)
%   width=1;
%   height=0.8;
%   subplot('position',[0. 0.14 width height])
%   m_proj('mercator',...
%          'lon',[domaxis(1) domaxis(2)],...
%          'lat',[domaxis(3) domaxis(4)]);
%   m_pcolor(Lonr,Latr,h.*themask);
%   shading flat
%   caxis(colaxis)
%   hold on
%   [C1,h1]=m_contour(Lonr,Latr,h,[hmin 100 200 500 1000 2000 4000],'k');
%   clabel(C1,h1,'LabelSpacing',1000,'Rotation',0,'Color','r')
%   if ~isempty(coastfileplot)
%     m_usercoast(coastfileplot,'color','r');
%     %m_usercoast(coastfileplot,'speckle','color','r');
%   else
%     m_gshhs_l('color','r');
%     m_gshhs_l('speckle','color','r');
%   end
%   m_grid('box','fancy',...
%          'xtick',5,'ytick',5,'tickdir','out',...
%          'fontsize',7);
%   hold off
% end
% warning on
% end
