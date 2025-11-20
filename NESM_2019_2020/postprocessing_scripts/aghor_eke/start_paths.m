clear all
close all
%------------------------------------
% on pace
main_dir_path='/storage/home/hcoda1/9/paghor3/';
aghor_extras_path='/storage/home/hcoda1/9/paghor3/aghor_extras/'; % change according to the system
%------------------------------------
% on prometheus
% main_dir_path='/home/aghor/';
% aghor_extras_path='/home/aghor/aghor/GT/GT_project/aghor_extras/'; % change according to the system
%------------------------------------
plots_path='../../plots/';
addpath([aghor_extras_path, 'export_fig/']);
path_str ='../../data/';
% path_str = strrep(path_str, '\','/');
%-------------------------------------------------------
% addpath to m_map library
% addpath('E:/GT_project/m_map1.4/m_map')
% addpath('/home/aghor/aghor/GT/GT_project/m_map') % on prometheus
addpath([aghor_extras_path, 'm_map/']); % on pace cluster
%------------------------------------
addpath([aghor_extras_path, 'cmap_manual/']); % manual colormaps with white in the middle
%------------------------------------
addpath([aghor_extras_path, 'gsw_matlab_v3_06_16/']); % on pace cluster
addpath([aghor_extras_path, 'gsw_matlab_v3_06_16/library']);
%------------------------------------
dirHR = [path_str, 'output/'];
grid_path = [path_str, 'edit_grid/'];
addpath(dirHR);
addpath(grid_path);
%------------------------------------
type = 'r';
coef = 1;
gridfile = 'NESM_grd.nc'; 
% s_rho = ncread(gridfile, 's_rho');
theta_s=3; % ncread(gridfile, 'theta_s');
theta_b=5; % ncread(gridfile, 'theta_b');
maskr=pagetranspose(ncread(gridfile, 'mask_rho'));
masku=pagetranspose(ncread(gridfile, 'mask_u'));
maskv=pagetranspose(ncread(gridfile, 'mask_v'));

Vtransform=2; % ncread(gridfile, 'Vtransform');
NumLayers = 100; 
lat_rho = pagetranspose(ncread(gridfile, 'lat_rho'));
lon_rho = pagetranspose(ncread(gridfile, 'lon_rho'));
depth = pagetranspose(ncread(gridfile, 'h'));
hc = 30; % ncread(gridfile, 'hc');

lon_rho_vec = squeeze(lon_rho(1, :));
lat_rho_vec = squeeze(lat_rho(:, 1));

pm = pagetranspose(ncread(gridfile, 'pm'));
pn = pagetranspose(ncread(gridfile, 'pn'));

mask_rho = pagetranspose(ncread(gridfile, 'mask_rho'));
mask_u = pagetranspose(ncread(gridfile, 'mask_u'));
mask_v = pagetranspose(ncread(gridfile, 'mask_v'));
%--------------------------------------------------------------------------
% taken from croco_tools/start.m
% main_dir_path = '/storage/home/hcoda1/9/paghor3';
tools_path = [main_dir_path, 'croco/croco_tools/'];
myutilpath=[tools_path,'UTILITIES/'];
addpath(tools_path);
addpath(myutilpath);
addpath([tools_path,'Diagnostic_tools'])
addpath([tools_path,'Preprocessing_tools'])
addpath([tools_path,'Visualization_tools'])
%-------------------------------------------------------
%-------------------------------------------------------
%
addpath([myutilpath,'mexcdf/mexnc'])   % 32 and 64 bits version of mexnc 
%
% - If these directories are already in your matlab native path, 
% you can comment these lines
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncsource'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/nctype'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncutility'])
%
% Use of built in opendap libraries (no loaddap) - S. Illig 2015 
%
  addpath([tools_path,'Opendap_tools_no_loaddap'])
%
%-------------------------------------------------------
%
% Use of loaddap  (older versions of matlab)
%
  addpath([tools_path,'Opendap_tools'])
%-------------------------------------------------------
%------------------------------------
%Extract NC file variables
time    = 0;
fname   = sprintf('NESM_2019_2020_avg.%05d.nc',time); hisfile = fname;
temp = ncread([dirHR fname],'temp');
% Calculate Depth
tindex = 1;
[zr]=get_depths(fname,gridfile,tindex,type);
[zw]=get_depths(fname,gridfile,tindex,'w');
[N, M, L] = size(zr);    
Ny = M; Nx = L;
% zr = shiftdim(zr,1);
%------------------------------------
% Atlantis II lat-lon lims
atl_lonmin = -63.5;   % Minimum longitude of Atlantis II  [degree east]
atl_lonmax = -62.84;   % Maximum longitude of Atlantis II [degree east]
atl_latmin = 38;   % Minimum latitude of Atlantis II [degree north]
atl_latmax = 38.9;   % Maximum latitude of Atlantis II [degree north]
%-------------------------------------------------------
rho0 = 1025; % ref density in kg/m^3
Omega = 7.2921e-5; % rotation rate of Earth s^-1
g = 9.81; % m/s^2
f0 = 0.909*1e-4; % Coriolis freq ;
% f0 = 2*Omega*sin(lat_rho_vec(Ny/2, 1)); % take average f0 over latitudes
N0 = 1e-3; % Brunt-Vaisala frequency (order of magnitude, s^-1)
%-------------------------------------------------------
% box limits for horizontal averages
lon_idx_min = 120;
lon_idx_max = 195;
lat_idx_min = 55;
lat_idx_max = 155;
%-------------------------------------------------------
