clear all
close all
%------------------------------------
% on pace
% main_dir_path='/storage/home/hcoda1/9/paghor3/';
% aghor_extras_path='/storage/home/hcoda1/9/paghor3/aghor_extras/'; % change according to the system
%------------------------------------
% on prometheus
% main_dir_path='/home/aghor/';
% aghor_extras_path='/home/aghor/aghor/GT/GT_project/aghor_extras/'; % change according to the system
%------------------------------------
% on geometer
main_dir_path='F:\prometheus_backup\Aghor\GT\GT_project\';
aghor_extras_path = 'F:\prometheus_backup\Aghor\GT\GT_project\aghor_extras\';
%------------------------------------
plots_path='./delta_fields_2d';
addpath([aghor_extras_path, 'export_fig/']);
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
addpath([tools_path,'Diagnostic_tools/Transport'])
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
nesm_1km_gridfile = '../../../NESM_2019_2020/data/edit_grid/NESM_grd.nc';
nesm_5km_gridfile = '../../../NESM_2019_2020_5km/data/edit_grid/NESM_grd_5km.nc';
nesm_5km_to_1km_gridfile = '../../../NESM_2019_2020_5km_to_1km/data/edit_grid/NESM_grd_5km_to_1km.nc';
sm3_1km_gridfile = '../../../SM3_2019_2020/data/edit_grid/SM3_grd.nc';
sm1_1km_gridfile = '../../../SM1_2019_2020/data/edit_grid/SM1_grd.nc';
%------------------------------------
theta_s=3; % ncread(gridfile, 'theta_s');
theta_b=5; % ncread(gridfile, 'theta_b');
Vtransform=2; % ncread(gridfile, 'Vtransform');
NumLayers = 100; 
%------------------------------------
nesm_1km_lat_rho = pagetranspose(ncread(nesm_1km_gridfile, 'lat_rho'));
nesm_1km_lon_rho = pagetranspose(ncread(nesm_1km_gridfile, 'lon_rho'));

nesm_5km_lat_rho = pagetranspose(ncread(nesm_5km_gridfile, 'lat_rho'));
nesm_5km_lon_rho = pagetranspose(ncread(nesm_5km_gridfile, 'lon_rho'));

nesm_5km_to_1km_lat_rho = pagetranspose(ncread(nesm_5km_to_1km_gridfile, 'lat_rho'));
nesm_5km_to_1km_lon_rho = pagetranspose(ncread(nesm_5km_to_1km_gridfile, 'lon_rho'));

sm3_1km_lat_rho = pagetranspose(ncread(sm3_1km_gridfile, 'lat_rho'));
sm3_1km_lon_rho = pagetranspose(ncread(sm3_1km_gridfile, 'lon_rho'));

sm1_1km_lat_rho = pagetranspose(ncread(sm1_1km_gridfile, 'lat_rho'));
sm1_1km_lon_rho = pagetranspose(ncread(sm1_1km_gridfile, 'lon_rho'));

% for plotting
lon_rho_1km = nesm_1km_lon_rho;
lat_rho_1km = nesm_1km_lat_rho;

lon_rho_5km = nesm_5km_lon_rho;
lat_rho_5km = nesm_5km_lat_rho;

lon_rho_5km_to_1km = nesm_5km_to_1km_lon_rho;
lat_rho_5km_to_1km = nesm_5km_to_1km_lat_rho;

lon_rho_vec_1km = squeeze(nesm_1km_lon_rho(1, :));
lat_rho_vec_1km = squeeze(nesm_1km_lat_rho(:, 1));

lon_rho_vec_5km = squeeze(nesm_5km_lon_rho(1, :));
lat_rho_vec_5km = squeeze(nesm_5km_lat_rho(:, 1));

lon_rho_vec_5km_to_1km = squeeze(nesm_5km_to_1km_lon_rho(1, :));
lat_rho_vec_5km_to_1km = squeeze(nesm_5km_to_1km_lat_rho(:, 1));

nesm_1km_depth = pagetranspose(ncread(nesm_1km_gridfile, 'h'));
nesm_5km_depth = pagetranspose(ncread(nesm_5km_gridfile, 'h'));
nesm_5km_to_1km_depth = pagetranspose(ncread(nesm_5km_to_1km_gridfile, 'h'));

sm3_1km_depth = pagetranspose(ncread(sm3_1km_gridfile, 'h'));
sm1_1km_depth = pagetranspose(ncread(sm1_1km_gridfile, 'h'));

%------------------------------------
f0 = 0.909e-4; % Coriolis freq.
%------------------------------------