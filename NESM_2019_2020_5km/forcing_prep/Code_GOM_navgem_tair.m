%%%%%%% This code is for large forcing file generation %%%%%%%%%
tnum_ori = datenum('1900-12-31 00:00:00');
t_min = tnum_start-tnum_ori;
t_max = tnum_end-tnum_ori;
var_in = 'airtmp';

varname = 'tair';
longname_var = 'surface air temperature';
unit_var = 'Celsius';

disp(['Reading data ...'])

time1 = ncread(tairfile1,'MT');
II = find(time1>=t_min & time1<=t_max);
II = II(1:dt:end);
time11 = time1(II);
Nt1 = length(II);
tt1 = ncread(tairfile1,var_in,[1 1 II(1)],[Inf Inf length(II)],[1 1 dt]);

if time11(end)<t_max
  time2 = ncread(tairfile2,'MT');
  II = find(time2>=(time11(end)+dt/24) & time2<=t_max);
  II = II(1:dt:end);
  time2 = time2(II);
  Nt2 = length(II);
  tt2 = ncread(tairfile2,var_in,[1 1 II(1)],[Inf Inf length(II)]);
else
  Nt2 = 0;
end
time1 = time11;

NNt = size(tt1);
frc_time = nan(Nt1+Nt2,1);
tt = nan(NNt(1),NNt(2),Nt1+Nt2);

frc_time(1:Nt1) = time1;
tt(:,:,1:Nt1) = tt1;
clear tt1
if Nt2 > 0
  frc_time(Nt1+1:end) = time2;
  tt(:,:,Nt1+1:end) = tt2;
  clear tt2
end
frc_time = frc_time + tnum_ori - datenum(Yorig,1,1);
tt = tt - 273.15; %K->degC
%frc_time = frc_time(1:dt:end);
%tt = tt(:,:,1:dt:end);

if ~isequal(frc_time,bulk_time)
  disp(['ERROR: time of ',varname,' does not match bulk_time!!!'])
  return
end

lon0 = ncread(tairfile1,'Longitude')-360;
lat0 = ncread(tairfile1,'Latitude');
[lat0,lon0]=meshgrid(lat0,lon0);

tti = nan(Lp,Mp,length(frc_time));
disp(['Interpolation ...'])
for it = 1:(length(frc_time))
  %disp(['  Day: ',num2str(frc_time(it))])
  tti(:,:,it) = interp2(lat0,lon0,squeeze(tt(:,:,it)),latr,lonr);
end

disp(['Saving data ...'])
ncid = netcdf.open(blkout,'WRITE');
netcdf.reDef(ncid);

varid2=netcdf.defVar(ncid,varname,'NC_FLOAT',[xi_rho eta_rho blkt]);

netcdf.putAtt(ncid,varid2,'long_name',longname_var)
netcdf.putAtt(ncid,varid2,'units',unit_var)
netcdf.putAtt(ncid,varid2,'positive',pos_var)

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid2,tti);

netcdf.close(ncid)

tair = tti; % save tair for the calculation of rhum

clear tti
clear tt
