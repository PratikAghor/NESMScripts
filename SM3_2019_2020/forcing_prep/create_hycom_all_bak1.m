function create_hycom_bry_all(in_begin,in_end,fname_out,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the OGCM file from hycom output that is ready to make crocco
% bry files
% by sdxmonkey
% 2022-12-05
%
%  modified from create_SODA
%  Further Information:  
%  http://www.croco-ocean.org
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% load informations from orignal hycom output file
file_in = [in_begin,'water_temp',in_end];
lonT = ncread(file_in,'lon');
latT = ncread(file_in,'lat');
depth = -ncread(file_in,'depth');
time = ncread(file_in,'time') + datenum('2000-01-01 00:00:00') - datenum(Yorig,1,1);
temp = ncread(file_in,'water_temp');
max(temp(:))
min(temp(:))
file_in = [in_begin,'salinity',in_end];
time_tmp = ncread(file_in,'time') + datenum('2000-01-01 00:00:00') - datenum(Yorig,1,1);
if ~isequal(time_tmp,time)
    disp(['ERROR: time of salinity is different from temperature!'])
    return
end
salt = ncread(file_in,'salinity');
max(salt(:))
min(salt(:))
file_in = [in_begin,'water_u',in_end];
time_tmp = ncread(file_in,'time') + datenum('2000-01-01 00:00:00') - datenum(Yorig,1,1);
if ~isequal(time_tmp,time)
    disp(['ERROR: time of U is different from temperature!'])
    return
end
lonU = ncread(file_in,'lon');
latU = ncread(file_in,'lat');
u = ncread(file_in,'water_u');
max(u(:))
min(u(:))
file_in = [in_begin,'water_v',in_end];
time_tmp = ncread(file_in,'time') + datenum('2000-01-01 00:00:00') - datenum(Yorig,1,1);
if ~isequal(time_tmp,time)
    disp(['ERROR: time of V is different from temperature!'])
    return
end
lonV = ncread(file_in,'lon');
latV = ncread(file_in,'lat');
v = ncread(file_in,'water_v');
max(v(:))
min(v(:))
file_in = [in_begin,'surf_el',in_end];
time_tmp = ncread(file_in,'time') + datenum('2000-01-01 00:00:00') - datenum(Yorig,1,1);
if ~isequal(time_tmp,time)
    disp(['ERROR: time of SSH is different from temperature!'])
    return
end
ssh = ncread(file_in,'surf_el');
max(ssh(:))
min(ssh(:))
% temp(temp > 1e30) = NaN;
% salt(salt > 1e30) = NaN;
% u(u > 1e30) = NaN;
% v(v > 1e30) = NaN;
% ssh(ssh>1e30) = NaN;

%missval=NaN;
disp('    Create the OGCM file')
ncid=netcdf.create(fname_out,'64bit_offset');
%redef(nc);
dimlonT=netcdf.defDim(ncid,'lonT',length(lonT));
dimlatT=netcdf.defDim(ncid,'latT',length(latT));
dimlonU=netcdf.defDim(ncid,'lonU',length(lonU));
dimlatU=netcdf.defDim(ncid,'latU',length(latU));
dimlonV=netcdf.defDim(ncid,'lonV',length(lonV));
dimlatV=netcdf.defDim(ncid,'latV',length(latV));
dimdepth=netcdf.defDim(ncid,'depth',length(depth));
dimtime=netcdf.defDim(ncid,'time',length(time));

idtemp = netcdf.defVar(ncid,'temp','float',[dimlonT dimlatT dimdepth dimtime]);
netcdf.putAtt(ncid,idtemp,'long_name','TEMPERATURE');
netcdf.putAtt(ncid,idtemp,'units','deg. C');

idsalt = netcdf.defVar(ncid,'salt','float',[dimlonT dimlatT dimdepth dimtime]);
netcdf.putAtt(ncid,idsalt,'long_name','SALINITY');
netcdf.putAtt(ncid,idsalt,'units','ppt');

idu = netcdf.defVar(ncid,'u','float',[dimlonU dimlatU dimdepth dimtime]);
netcdf.putAtt(ncid,idu,'long_name','ZONAL VELOCITY');
netcdf.putAtt(ncid,idu,'units','m/sec');

idv = netcdf.defVar(ncid,'v','float',[dimlonV dimlatV dimdepth dimtime]);
netcdf.putAtt(ncid,idv,'long_name','MERIDIONAL VELOCITY');
netcdf.putAtt(ncid,idv,'units','m/sec');

idubar = netcdf.defVar(ncid,'ubar','float',[dimlonU dimlatU dimtime]);
netcdf.putAtt(ncid,idubar,'long_name','ZONAL BAROTROPIC VELOCITY');
netcdf.putAtt(ncid,idubar,'units','m/sec');

idvbar = netcdf.defVar(ncid,'vbar','float',[dimlonV dimlatV dimtime]);
netcdf.putAtt(ncid,idvbar,'long_name','MERIDIONAL BAROTROPIC VELOCITY');
netcdf.putAtt(ncid,idvbar,'units','m/sec');

idssh = netcdf.defVar(ncid,'ssh','float',[dimlonT dimlatT dimtime]);
netcdf.putAtt(ncid,idssh,'long_name','SEA LEVEL HEIGHT');
netcdf.putAtt(ncid,idssh,'units','m');

idlonT = netcdf.defVar(ncid,'lonT','double',[dimlonT]);
netcdf.putAtt(ncid,idlonT,'units','degrees_east');

idlatT = netcdf.defVar(ncid,'latT','double',[dimlatT]);
netcdf.putAtt(ncid,idlatT,'units','degrees_north');


idlonU = netcdf.defVar(ncid,'lonU','double',[dimlonU]);
netcdf.putAtt(ncid,idlonU,'units','degrees_east');

idlatU = netcdf.defVar(ncid,'latU','double',[dimlatU]);
netcdf.putAtt(ncid,idlatU,'units','degrees_north');

idlonV = netcdf.defVar(ncid,'lonV','double',[dimlonV]);
netcdf.putAtt(ncid,idlonV,'units','degrees_east');

idlatV = netcdf.defVar(ncid,'latV','double',[dimlatV]);
netcdf.putAtt(ncid,idlatV,'units','degrees_north');

iddepth = netcdf.defVar(ncid,'depth','double',[dimdepth]);
netcdf.putAtt(ncid,iddepth,'units','meters');

idtime = netcdf.defVar(ncid,'time','double',[dimtime]);
netcdf.putAtt(ncid,idtime,'units',['days since 1-Jan-',num2str(Yorig),' 00:00:0.0']);
netcdf.endDef(ncid);
netcdf.close(ncid)
%
% File the file
%
disp('    Fill the OGCM file')
nc=netcdf(fname_out,'write')
nc{'depth'}(:)=depth;
nc{'latT'}(:)=latT;
nc{'lonT'}(:)=lonT;
nc{'latU'}(:)=latU;
nc{'lonU'}(:)=lonU;
nc{'latV'}(:)=latV;
nc{'lonV'}(:)=lonV;
%
for tndx=1:length(time)
%
nc{'time'}(tndx)=time(tndx);
%
if length(time)==1
  nc{'ssh'}(tndx,:,:)=ssh;
  u1=u;
  v1=v;
  nc{'u'}(tndx,:,:,:)=u1;
  nc{'v'}(tndx,:,:,:)=v1;
  nc{'temp'}(tndx,:,:,:)=temp;
  nc{'salt'}(tndx,:,:,:)=salt;
else
  nc{'ssh'}(tndx,:,:)=squeeze(ssh(tndx,:,:));
  u1=squeeze(u(tndx,:,:,:));
  v1=squeeze(v(tndx,:,:,:));
  nc{'u'}(tndx,:,:,:)=u1;
  nc{'v'}(tndx,:,:,:)=v1;
  nc{'temp'}(tndx,:,:,:)=squeeze(temp(tndx,:,:,:));
  nc{'salt'}(tndx,:,:,:)=squeeze(salt(tndx,:,:,:));
end
%
% Compute the barotropic velocities
%
masku=isfinite(u1);
maskv=isfinite(v1);
u1(isnan(u1))=0;
v1(isnan(v1))=0;
dz=gradient(depth);
NZ=length(depth);
du=0*squeeze(u1(1,:,:));
zu=du;
dv=0*squeeze(v1(1,:,:));
zv=dv;
for k=1:NZ
  du=du+dz(k)*squeeze(u1(k,:,:));
  zu=zu+dz(k)*squeeze(masku(k,:,:));
  dv=dv+dz(k)*squeeze(v1(k,:,:));
  zv=zv+dz(k)*squeeze(maskv(k,:,:));
end
du(zu==0)=NaN;
dv(zv==0)=NaN;
zu(zu==0)=NaN;
zv(zv==0)=NaN;
ubar=du./zu;
vbar=dv./zv;
%
nc{'ubar'}(tndx,:,:)=ubar;
nc{'vbar'}(tndx,:,:)=vbar;
%
end
%
close(nc)
%
return
