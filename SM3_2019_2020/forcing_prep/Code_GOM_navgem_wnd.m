%%%%%%% This code is for large forcing file generation %%%%%%%%%
tnum_ori = datenum('1900-12-31 00:00:00');
t_min = tnum_start-tnum_ori;
t_max = tnum_end-tnum_ori;
disp(['Reading data ...'])

time1 = ncread(wndfile1,'MT');
II = find(time1>=t_min & time1<=t_max);
II = II(1:dt:end);
time11 = time1(II);
Nt1 = length(II)
uu1 = single(ncread(wndfile1,'wndewd',[1 1 II(1)],[Inf Inf length(II)],[1 1 dt]));
vv1 = single(ncread(wndfile1,'wndnwd',[1 1 II(1)],[Inf Inf length(II)],[1 1 dt])); 

if time11(end)<t_max
  time2 = ncread(wndfile2,'MT');
  II = find(time2 > (time11(end)) & time2<=t_max);
  II = II(1:dt:end);
  time2 = time2(II);
  Nt2 = length(II)
  uu2 = single(ncread(wndfile2,'wndewd',[1 1 II(1)],[Inf Inf length(II)],[1 1 dt]));
  vv2 = single(ncread(wndfile2,'wndnwd',[1 1 II(1)],[Inf Inf length(II)],[1 1 dt]));
else
  Nt2 = 0;
end
time1 = time11;

NNu = size(uu1);
NNv = size(vv1);
bulk_time = nan(Nt1+Nt2,1);
uu = nan(NNu(1),NNu(2),Nt1+Nt2);
vv = nan(NNv(1),NNv(2),Nt1+Nt2);

bulk_time(1:Nt1) = time1;
uu(:,:,1:Nt1) = double(uu1);
vv(:,:,1:Nt1) = double(vv1);
clear uu1 vv1
if Nt2 > 0
  bulk_time(Nt1+1:end) = time2;
  uu(:,:,Nt1+1:end) = double(uu2);
  vv(:,:,Nt1+1:end) = double(vv2);
  clear uu2 vv2
end
bulk_time = bulk_time + tnum_ori - datenum(Yorig,1,1);
% hourly data -> 3-hour data
%bulk_time = bulk_time(1:dt:end);
%uu = uu(:,:,1:dt:end);
%vv = vv(:,:,1:dt:end);

lon0 = ncread(wndfile1,'Longitude')-360;
lat0 = ncread(wndfile1,'Latitude');
[lat0,lon0]=meshgrid(lat0,lon0);

uwnd = nan(Lp,Mp,length(bulk_time));
vwnd = nan(Lp,Mp,length(bulk_time));
disp(['Interpolation ...'])
for it = 1:length(bulk_time)%(Nt1+Nt2)
  %disp(['  Day: ',num2str(bulk_time(it))])
  uwnd(:,:,it) = interp2(lat0,lon0,squeeze(uu(:,:,it)),latr,lonr);
  vwnd(:,:,it) = interp2(lat0,lon0,squeeze(vv(:,:,it)),latr,lonr);
end
wspd=sqrt(uwnd.^2+vwnd.^2);

disp(['Saving data ...'])
ncid = netcdf.open(blkout,'WRITE');
netcdf.reDef(ncid);
blkt=netcdf.defDim(ncid,'bulk_time',length(bulk_time));

varid1=netcdf.defVar(ncid,'bulk_time','NC_FLOAT',[blkt]);
varid2=netcdf.defVar(ncid,'uwnd','NC_FLOAT',[xi_rho eta_rho blkt]);
varid3=netcdf.defVar(ncid,'vwnd','NC_FLOAT',[xi_rho eta_rho blkt]);
varid4=netcdf.defVar(ncid,'wspd','NC_FLOAT',[xi_rho eta_rho blkt]);

netcdf.putAtt(ncid,varid1,'long_name','bulk formulation execution time')
netcdf.putAtt(ncid,varid1,'units','days')

netcdf.putAtt(ncid,varid2,'long_name','u-wind')
netcdf.putAtt(ncid,varid2,'units','m/s')

netcdf.putAtt(ncid,varid3,'long_name','u-wind')
netcdf.putAtt(ncid,varid3,'units','m/s')

netcdf.putAtt(ncid,varid4,'long_name','wind speed 10m')
netcdf.putAtt(ncid,varid4,'units','m s-1')

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid1,bulk_time);
netcdf.putVar(ncid,varid2,uwnd);
netcdf.putVar(ncid,varid3,vwnd);
netcdf.putVar(ncid,varid4,wspd);

netcdf.close(ncid)

clear uwnd
clear vwnd
clear uu
clear vv
clear wspd
