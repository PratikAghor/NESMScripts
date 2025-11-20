function create_bryfile(bryname,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,cycle,clobber,vtransform);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer 
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')
if nargin < 12
    disp([' NO VTRANSFORM parameter found'])
    disp([' USE TRANSFORM default value vtransform = 1'])
    vtransform = 1; 
end
disp([' VTRANSFORM = ',num2str(vtransform)])
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
close(nc);
hmin=min(min(h(maskr==1)));
if vtransform ==1;
  if hc > hmin
    error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
  end
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ; 
history = 'CROCO' ;
ncid=netcdf.create(bryname,'64BIT_OFFSET');
% aghoredit: close and open in reDef mode
% netcdf.close(ncid);
% netcdf.open(bryname, 'WRITE');
% netcdf.reDef(ncid);
%%result = redef(nc);
%
%  Create dimensions
%
xi_u=netcdf.defDim(ncid,'xi_u',L);
xi_v=netcdf.defDim(ncid,'xi_v',Lp);
xi_rho=netcdf.defDim(ncid,'xi_rho',Lp);
eta_u=netcdf.defDim(ncid,'eta_u',Mp);
eta_v=netcdf.defDim(ncid,'eta_v',M);
eta_rho=netcdf.defDim(ncid,'eta_rho',Mp);
s_rho=netcdf.defDim(ncid,'s_rho',N);
s_w=netcdf.defDim(ncid,'s_w',Np);
tracer=netcdf.defDim(ncid,'tracer',2);
bry_t=netcdf.defDim(ncid,'bry_time',length(time));
tclm_t=netcdf.defDim(ncid,'tclm_time',length(time));
temp_t=netcdf.defDim(ncid,'temp_time',length(time));
sclm_t=netcdf.defDim(ncid,'sclm_time',length(time));
salt_t=netcdf.defDim(ncid,'salt_time',length(time));
uclm_t=netcdf.defDim(ncid,'uclm_time',length(time));
vclm_t=netcdf.defDim(ncid,'vclm_time',length(time));
v2d_t=netcdf.defDim(ncid,'v2d_time',length(time));
v3d_t=netcdf.defDim(ncid,'v3d_time',length(time));
ssh_t=netcdf.defDim(ncid,'ssh_time',length(time));
zeta_t=netcdf.defDim(ncid,'zeta_time',length(time));
one=netcdf.defDim(ncid,'one',1);
%
%  Create variables and attributes
%
spherical = netcdf.defVar(ncid,'spherical','char',[one]);
netcdf.putAtt(ncid,spherical,'long_name','grid type logical switch');
netcdf.putAtt(ncid,spherical,'flag_values','T, F');
netcdf.putAtt(ncid,spherical,'flag_meanings','spherical Cartesian');

%
Vtransform = netcdf.defVar(ncid,'Vtransform','int',[one]);
netcdf.putAtt(ncid,Vtransform,'long_name','vertical terrain-following transformation equation');
%
Vstretching = netcdf.defVar(ncid,'Vstretching','int',[one]);
netcdf.putAtt(ncid,Vstretching,'long_name','vertical terrain-following stretching function');
%
tstart = netcdf.defVar(ncid,'tstart','double',[one]);
netcdf.putAtt(ncid,tstart,'long_name','start processing day');
netcdf.putAtt(ncid,tstart,'units','day');
%
tend = netcdf.defVar(ncid,'tend','double',[one]);
netcdf.putAtt(ncid,tend,'long_name','end processing day');
netcdf.putAtt(ncid,tend,'units','day');
%
theta_s_ind = netcdf.defVar(ncid,'theta_s','double',[one]);
netcdf.putAtt(ncid,theta_s_ind,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(ncid,theta_s_ind,'units','nondimensional');
%
theta_b_ind = netcdf.defVar(ncid,'theta_b','double',[one]);
netcdf.putAtt(ncid,theta_b_ind,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(ncid,theta_b_ind,'units','nondimensional');
%
Tcline_ind = netcdf.defVar(ncid,'Tcline','double',[one]);
netcdf.putAtt(ncid,Tcline_ind,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(ncid,Tcline_ind,'units','meter');
%
hc_ind = netcdf.defVar(ncid,'hc','double',[one]);
netcdf.putAtt(ncid,hc_ind,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(ncid,hc_ind,'units','meter');
%
sc_r_ind = netcdf.defVar(ncid,'sc_r','double',[s_rho]);
netcdf.putAtt(ncid,sc_r_ind,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(ncid,sc_r_ind,'valid_min',-1.);
netcdf.putAtt(ncid,sc_r_ind,'valid_min',0.);
netcdf.putAtt(ncid,sc_r_ind,'positive','up');
if (vtransform == 1)
    netcdf.putAtt(ncid,sc_r_ind,'standard_name','ocena_s_coordinate_g1');
elseif (vtransform == 2)
    netcdf.putAtt(ncid,sc_r_ind,'standard_name','ocena_s_coordinate_g2');
end
netcdf.putAtt(ncid,sc_r_ind,'formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
%
sc_w_ind = netcdf.defVar(ncid,'sc_w','double',[s_w]);
netcdf.putAtt(ncid,sc_w_ind,'long_name','S-coordinate at W-points');
netcdf.putAtt(ncid,sc_w_ind,'valid_min',-1.);
netcdf.putAtt(ncid,sc_w_ind,'valid_min',0.);
netcdf.putAtt(ncid,sc_w_ind,'positive','up');
if (vtransform == 1)
    netcdf.putAtt(ncid,sc_w_ind,'standard_name','ocena_s_coordinate_g1');
elseif (vtransform == 2)
    netcdf.putAtt(ncid,sc_w_ind,'standard_name','ocena_s_coordinate_g2');
end
netcdf.putAtt(ncid,sc_w_ind,'formula_terms','s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc');
%
Cs_r_ind = netcdf.defVar(ncid,'Cs_r','double',[s_rho]);
netcdf.putAtt(ncid,Cs_r_ind,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(ncid,Cs_r_ind,'units','nondimensional');
netcdf.putAtt(ncid,Cs_r_ind,'valid_min',-1.);
netcdf.putAtt(ncid,Cs_r_ind,'valid_min',0.);
%
Cs_w_ind = netcdf.defVar(ncid,'Cs_w','double',[s_w]);
netcdf.putAtt(ncid,Cs_w_ind,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(ncid,Cs_w_ind,'units','nondimensional');
netcdf.putAtt(ncid,Cs_w_ind,'valid_min',-1.);
netcdf.putAtt(ncid,Cs_w_ind,'valid_min',0.);
%
bry_time = netcdf.defVar(ncid,'bry_time','double',[bry_t]);
netcdf.putAtt(ncid,bry_time,'long_name','time for boundary climatology');
netcdf.putAtt(ncid,bry_time,'units','day');
netcdf.putAtt(ncid,bry_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,bry_time,'cycle_length',cycle);
%
tclm_time = netcdf.defVar(ncid,'tclm_time','double',[tclm_t]);
netcdf.putAtt(ncid,tclm_time,'long_name','time for temperature climatology');
netcdf.putAtt(ncid,tclm_time,'units','day');
netcdf.putAtt(ncid,tclm_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,tclm_time,'cycle_length',cycle);
%
temp_time = netcdf.defVar(ncid,'temp_time','double',[temp_t]);
netcdf.putAtt(ncid,temp_time,'long_name','time for temperature climatology');
netcdf.putAtt(ncid,temp_time,'units','day');
netcdf.putAtt(ncid,temp_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,temp_time,'cycle_length',cycle);
%
sclm_time = netcdf.defVar(ncid,'sclm_time','double',[sclm_t]);
netcdf.putAtt(ncid,sclm_time,'long_name','time for salinity climatology');
netcdf.putAtt(ncid,sclm_time,'units','day');
netcdf.putAtt(ncid,sclm_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,sclm_time,'cycle_length',cycle);
%
salt_time = netcdf.defVar(ncid,'salt_time','double',[salt_t]);
netcdf.putAtt(ncid,salt_time,'long_name','time for salinity climatology');
netcdf.putAtt(ncid,salt_time,'units','day');
netcdf.putAtt(ncid,salt_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,salt_time,'cycle_length',cycle);
%
uclm_time = netcdf.defVar(ncid,'uclm_time','double',[uclm_t]);
netcdf.putAtt(ncid,uclm_time,'long_name','time climatological u');
netcdf.putAtt(ncid,uclm_time,'units','day');
netcdf.putAtt(ncid,uclm_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,uclm_time,'cycle_length',cycle);
%
vclm_time = netcdf.defVar(ncid,'vclm_time','double',[vclm_t]);
netcdf.putAtt(ncid,vclm_time,'long_name','time climatological v');
netcdf.putAtt(ncid,vclm_time,'units','day');
netcdf.putAtt(ncid,vclm_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,vclm_time,'cycle_length',cycle);
%
v2d_time = netcdf.defVar(ncid,'v2d_time','double',[v2d_t]);
netcdf.putAtt(ncid,v2d_time,'long_name','time for 2D velocity climatology');
netcdf.putAtt(ncid,v2d_time,'units','day');
netcdf.putAtt(ncid,v2d_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,v2d_time,'cycle_length',cycle);
%
v3d_time = netcdf.defVar(ncid,'v3d_time','double',[v3d_t]);
netcdf.putAtt(ncid,v3d_time,'long_name','time for 3D velocity climatology');
netcdf.putAtt(ncid,v3d_time,'units','day');
netcdf.putAtt(ncid,v3d_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,v3d_time,'cycle_length',cycle);
%
ssh_time = netcdf.defVar(ncid,'ssh_time','double',[ssh_t]);
netcdf.putAtt(ncid,ssh_time,'long_name','time for sea surface height');
netcdf.putAtt(ncid,ssh_time,'units','day');
netcdf.putAtt(ncid,ssh_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,ssh_time,'cycle_length',cycle);
%
zeta_time = netcdf.defVar(ncid,'zeta_time','double',[zeta_t]);
netcdf.putAtt(ncid,zeta_time,'long_name','time for sea surface height');
netcdf.putAtt(ncid,zeta_time,'units','day');
netcdf.putAtt(ncid,zeta_time,'calendar','360.0 days in every year');
netcdf.putAtt(ncid,zeta_time,'cycle_length',cycle);
%
if obc(1)==1
%
%   Southern boundary
%
  temp_south = netcdf.defVar(ncid,'temp_south','double',[xi_rho s_rho temp_t]);
  netcdf.putAtt(ncid,temp_south,'long_name','southern boundary potential temperature');
  netcdf.putAtt(ncid,temp_south,'units','Celsius');
  netcdf.putAtt(ncid,temp_south,'coordinates','lon_rho s_rho temp_time');
%
  salt_south = netcdf.defVar(ncid,'salt_south','double',[xi_rho s_rho salt_t]);
  netcdf.putAtt(ncid,salt_south,'long_name','southern boundary salinity');
  netcdf.putAtt(ncid,salt_south,'units','PSU');
  netcdf.putAtt(ncid,salt_south,'coordinates','lon_rho s_rho salt_time');
%
  u_south = netcdf.defVar(ncid,'u_south','double',[xi_u s_rho v3d_t]);
  netcdf.putAtt(ncid,u_south,'long_name','southern boundary u-momentum component');
  netcdf.putAtt(ncid,u_south,'units','meter second-1');
  netcdf.putAtt(ncid,u_south,'coordinates','lon_u s_rho uclm_time');
%
  v_south = netcdf.defVar(ncid,'v_south','double',[xi_rho s_rho v3d_t]);
  netcdf.putAtt(ncid,v_south,'long_name','southern boundary v-momentum component');
  netcdf.putAtt(ncid,v_south,'units','meter second-1');
  netcdf.putAtt(ncid,v_south,'coordinates','lon_v s_rho vclm_time');
%
  ubar_south = netcdf.defVar(ncid,'ubar_south','double',[xi_u v2d_t]);
  netcdf.putAtt(ncid,ubar_south,'long_name','southern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid,ubar_south,'units','meter second-1');
  netcdf.putAtt(ncid,ubar_south,'coordinates','lon_u uclm_time');
%
  vbar_south = netcdf.defVar(ncid,'vbar_south','double',[xi_rho v2d_t]);
  netcdf.putAtt(ncid,vbar_south,'long_name','southern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid,vbar_south,'units','meter second-1');
  netcdf.putAtt(ncid,vbar_south,'coordinates','lon_v vclm_time');
%
  zeta_south = netcdf.defVar(ncid,'zeta_south','double',[xi_rho zeta_t]);
  netcdf.putAtt(ncid,zeta_south,'long_name','southern boundary sea surface height');
  netcdf.putAtt(ncid,zeta_south,'units','meter');
  netcdf.putAtt(ncid,zeta_south,'coordinates','lon_rho zeta_time');
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  temp_east = netcdf.defVar(ncid,'temp_east','double',[eta_rho s_rho temp_t]);
  netcdf.putAtt(ncid,temp_east,'long_name','eastern boundary potential temperature');
  netcdf.putAtt(ncid,temp_east,'units','Celsius');
  netcdf.putAtt(ncid,temp_east,'coordinates','lat_rho s_rho temp_time');
%
  salt_east = netcdf.defVar(ncid,'salt_east','double',[eta_rho s_rho salt_t]);
  netcdf.putAtt(ncid,salt_east,'long_name','eastern boundary salinity');
  netcdf.putAtt(ncid,salt_east,'units','PSU');
  netcdf.putAtt(ncid,salt_east,'coordinates','lat_rho s_rho salt_time');
%
  u_east = netcdf.defVar(ncid,'u_east','double',[eta_rho s_rho v3d_t]);
  netcdf.putAtt(ncid,u_east,'long_name','eastern boundary u-momentum component');
  netcdf.putAtt(ncid,u_east,'units','meter second-1');
  netcdf.putAtt(ncid,u_east,'coordinates','lat_u s_rho uclm_time');
%
  v_east = netcdf.defVar(ncid,'v_east','double',[eta_v s_rho v3d_t]);
  netcdf.putAtt(ncid,v_east,'long_name','eastern boundary v-momentum component');
  netcdf.putAtt(ncid,v_east,'units','meter second-1');
  netcdf.putAtt(ncid,v_east,'coordinates','lat_v s_rho vclm_time');
%
  ubar_east = netcdf.defVar(ncid,'ubar_east','double',[eta_rho v2d_t]);
  netcdf.putAtt(ncid,ubar_east,'long_name','eastern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid,ubar_east,'units','meter second-1');
  netcdf.putAtt(ncid,ubar_east,'coordinates','lat_u uclm_time');
%
  vbar_east = netcdf.defVar(ncid,'vbar_east','double',[eta_v v2d_t]);
  netcdf.putAtt(ncid,vbar_east,'long_name','eastern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid,vbar_east,'units','meter second-1');
  netcdf.putAtt(ncid,vbar_east,'coordinates','lat_v vclm_time');
%
  zeta_east = netcdf.defVar(ncid,'zeta_east','double',[eta_rho zeta_t]);
  netcdf.putAtt(ncid,zeta_east,'long_name','eastern boundary sea surface height');
  netcdf.putAtt(ncid,zeta_east,'units','meter');
  netcdf.putAtt(ncid,zeta_east,'coordinates','lat_rho zeta_time');
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  temp_north = netcdf.defVar(ncid,'temp_north','double',[xi_rho s_rho temp_t]);
  netcdf.putAtt(ncid,temp_north,'long_name','northern boundary potential temperature');
  netcdf.putAtt(ncid,temp_north,'units','Celsius');
  netcdf.putAtt(ncid,temp_north,'coordinates','lon_rho s_rho temp_time');
%
  salt_north = netcdf.defVar(ncid,'salt_north','double',[xi_rho s_rho salt_t]);
  netcdf.putAtt(ncid,salt_north,'long_name','northern boundary salinity');
  netcdf.putAtt(ncid,salt_north,'units','PSU');
  netcdf.putAtt(ncid,salt_north,'coordinates','lon_rho s_rho salt_time');
%
  u_north = netcdf.defVar(ncid,'u_north','double',[xi_u s_rho v3d_t]);
  netcdf.putAtt(ncid,u_north,'long_name','northern boundary u-momentum component');
  netcdf.putAtt(ncid,u_north,'units','meter second-1');
  netcdf.putAtt(ncid,u_north,'coordinates','lon_u s_rho uclm_time');
%
  v_north = netcdf.defVar(ncid,'v_north','double',[xi_rho s_rho v3d_t]);
  netcdf.putAtt(ncid,v_north,'long_name','northern boundary v-momentum component');
  netcdf.putAtt(ncid,v_north,'units','meter second-1');
  netcdf.putAtt(ncid,v_north,'coordinates','lon_v s_rho vclm_time');
%
  ubar_north = netcdf.defVar(ncid,'ubar_north','double',[xi_u v2d_t]);
  netcdf.putAtt(ncid,ubar_north,'long_name','northern boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid,ubar_north,'units','meter second-1');
  netcdf.putAtt(ncid,ubar_north,'coordinates','lon_u uclm_time');
%
  vbar_north = netcdf.defVar(ncid,'vbar_north','double',[xi_rho v2d_t]);
  netcdf.putAtt(ncid,vbar_north,'long_name','northern boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid,vbar_north,'units','meter second-1');
  netcdf.putAtt(ncid,vbar_north,'coordinates','lon_v vclm_time');
%
  zeta_north = netcdf.defVar(ncid,'zeta_north','double',[xi_rho zeta_t]);
  netcdf.putAtt(ncid,zeta_north,'long_name','northern boundary sea surface height');
  netcdf.putAtt(ncid,zeta_north,'units','meter');
  netcdf.putAtt(ncid,zeta_north,'coordinates','lon_rho zeta_time');
%
end
%
if obc(4)==1
%
%   Western boundary
%
  temp_west = netcdf.defVar(ncid,'temp_west','double',[eta_rho s_rho temp_t]);
  netcdf.putAtt(ncid,temp_west,'long_name','western boundary potential temperature');
  netcdf.putAtt(ncid,temp_west,'units','Celsius');
  netcdf.putAtt(ncid,temp_west,'coordinates','lat_rho s_rho temp_time');
%
  salt_west = netcdf.defVar(ncid,'salt_west','double',[eta_rho s_rho salt_t]);
  netcdf.putAtt(ncid,salt_west,'long_name','western boundary salinity');
  netcdf.putAtt(ncid,salt_west,'units','PSU');
  netcdf.putAtt(ncid,salt_west,'coordinates','lat_rho s_rho salt_time');
%
  u_west = netcdf.defVar(ncid,'u_west','double',[eta_rho s_rho v3d_t]);
  netcdf.putAtt(ncid,u_west,'long_name','western boundary u-momentum component');
  netcdf.putAtt(ncid,u_west,'units','meter second-1');
  netcdf.putAtt(ncid,u_west,'coordinates','lat_u s_rho uclm_time');
%
  v_west = netcdf.defVar(ncid,'v_west','double',[eta_v s_rho v3d_t]);
  netcdf.putAtt(ncid,v_west,'long_name','western boundary v-momentum component');
  netcdf.putAtt(ncid,v_west,'units','meter second-1');
  netcdf.putAtt(ncid,v_west,'coordinates','lat_v s_rho vclm_time');
%
  ubar_west = netcdf.defVar(ncid,'ubar_west','double',[eta_rho v2d_t]);
  netcdf.putAtt(ncid,ubar_west,'long_name','western boundary vertically integrated u-momentum component');
  netcdf.putAtt(ncid,ubar_west,'units','meter second-1');
  netcdf.putAtt(ncid,ubar_west,'coordinates','lat_u uclm_time');
%
  vbar_west = netcdf.defVar(ncid,'vbar_west','double',[eta_v v2d_t]);
  netcdf.putAtt(ncid,vbar_west,'long_name','western boundary vertically integrated v-momentum component');
  netcdf.putAtt(ncid,vbar_west,'units','meter second-1');
  netcdf.putAtt(ncid,vbar_west,'coordinates','lat_v vclm_time');
%
  zeta_west = netcdf.defVar(ncid,'zeta_west','double',[eta_rho zeta_t]);
  netcdf.putAtt(ncid,zeta_west,'long_name','western boundary sea surface height');
  netcdf.putAtt(ncid,zeta_west,'units','meter');
  netcdf.putAtt(ncid,zeta_west,'coordinates','lat_rho zeta_time');
%
end
%
%
% Create global attributes
%
v_glob = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,v_glob,'title',title);
netcdf.putAtt(ncid,v_glob,'date',date);
netcdf.putAtt(ncid,v_glob,'clim_file',bryname);
netcdf.putAtt(ncid,v_glob,'grd_file',grdname);
netcdf.putAtt(ncid,v_glob,'type',type);
netcdf.putAtt(ncid,v_glob,'history',history);
netcdf.endDef(ncid);
netcdf.close(ncid)
%
% Leave define mode
%
%%result = endef(nc);
%
% Compute S coordinates
%
[sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform);
%disp(['vtransform=',num2str(vtransform)])
%
% Write variables
%
nc = netcdf(bryname,'write');
nc{'spherical'}(:)='T';
nc{'Vtransform'}(:)=vtransform;
nc{'Vstretching'}(:)=1;
nc{'tstart'}(:) =  min([min(time) min(time) min(time)]); 
nc{'tend'}(:) =  max([max(time) max(time) max(time)]); 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'sc_r'}(:) = sc_r;
nc{'sc_w'}(:) = sc_w;
nc{'Cs_r'}(:) = Cs_r ; 
nc{'Cs_w'}(:) = Cs_w;
nc{'tclm_time'}(:) =  time; 
nc{'temp_time'}(:) =  time; 
nc{'sclm_time'}(:) =  time; 
nc{'salt_time'}(:) =  time; 
nc{'uclm_time'}(:) =  time; 
nc{'vclm_time'}(:) =  time; 
nc{'v2d_time'}(:) =   time; 
nc{'v3d_time'}(:) =   time; 
nc{'ssh_time'}(:) =   time;
nc{'zeta_time'}(:) =  time;
nc{'bry_time'}(:) =  time; 
if obc(1)==1
  nc{'u_south'}(:) =  0; 
  nc{'v_south'}(:) =  0; 
  nc{'ubar_south'}(:) =  0; 
  nc{'vbar_south'}(:) =  0; 
  nc{'zeta_south'}(:) =  0; 
  nc{'temp_south'}(:) =  0; 
  nc{'salt_south'}(:) =  0;
end 
if obc(2)==1
  nc{'u_east'}(:) =  0; 
  nc{'v_east'}(:) =  0; 
  nc{'ubar_east'}(:) =  0; 
  nc{'vbar_east'}(:) =  0; 
  nc{'zeta_east'}(:) =  0; 
  nc{'temp_east'}(:) =  0; 
  nc{'salt_east'}(:) =  0;
end 
if obc(3)==1
  nc{'u_north'}(:) =  0; 
  nc{'v_north'}(:) =  0; 
  nc{'ubar_north'}(:) =  0; 
  nc{'vbar_north'}(:) =  0; 
  nc{'zeta_north'}(:) =  0; 
  nc{'temp_north'}(:) =  0; 
  nc{'salt_north'}(:) =  0;
end 
if obc(4)==1
  nc{'u_west'}(:) =  0; 
  nc{'v_west'}(:) =  0; 
  nc{'ubar_west'}(:) =  0; 
  nc{'vbar_west'}(:) =  0; 
  nc{'zeta_west'}(:) =  0; 
  nc{'temp_west'}(:) =  0; 
  nc{'salt_west'}(:) =  0;
end 
close(nc)
return


