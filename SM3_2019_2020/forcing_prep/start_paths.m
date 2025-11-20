clear all
close all
%------------------------------------
% on pace
% main_dir_path='/storage/home/hcoda1/9/paghor3/';
% aghor_extras_path='/storage/home/hcoda1/9/paghor3/aghor_extras/'; % change according to the system
%------------------------------------
% on prometheus
main_dir_path='/home/aghor/';
aghor_extras_path='/home/aghor/aghor/GT/GT_project/aghor_extras/'; % change according to the system
% on geometer
% main_dir_path='F:\prometheus_backup\Aghor\GT\GT_project\';
% aghor_extras_path = 'F:\prometheus_backup\Aghor\GT\GT_project\aghor_extras\';
%------------------------------------
% plots_path='../../plots/';
addpath([aghor_extras_path, 'export_fig/']);
path_str ='../data/';
%-------------------------------------------------------
% addpath to m_map library
% addpath('E:/GT_project/m_map1.4/m_map')
% addpath('/home/aghor/aghor/GT/GT_project/m_map') % on prometheus
addpath([aghor_extras_path, 'm_map/']); % on pace cluster
%------------------------------------
addpath([aghor_extras_path, 'cmap_manual/']); % manual colormaps with white in the middle
%------------------------------------
addpath([aghor_extras_path, 'gsw_matlab_v3_06_16/']); % on pace cluster
%------------------------------------
dirHR = [path_str, 'output/'];
grid_path = [path_str, 'edit_grid/'];
addpath(dirHR);
addpath(grid_path);
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
%-------------------------------------------------------
%-------------------------------------------------------
