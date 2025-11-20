% save annual mean of rho, w to calculate VEBF about this mean
%-----------------------------
clear all
close all
start_paths
%%
Mu = M; Lu = L-1; Mv = M-1; Lv = L;
% read files and calculate mean (bar) values of u, v
indxRange = 952:3877; % What time indices do you need?
nt0=indxRange(1);
[~, Nt] = size(indxRange);
time_arr = zeros(Nt, 1);
% ub=zeros(N, Mu, Lu);
% vb=zeros(N, Mv, Lv);
wb=zeros(N, M, L);
rhob=zeros(N, M, L);
% zetab=zeros(M, L);

% size(ub)
% size(u2rho_3d(ub))

for nt = 1:Nt
    nt
    % sprintf('nt = %d', nt);
    fname = sprintf([dirHR, 'NESM_2019_2020_5km_avg.%05d.nc'], indxRange(1, nt));
    hisfile = fname;
    nc=netcdf(fname);
    time = ncread(fname, 'time');
    time_arr(nt) = time;
    
    rhotmp = zeros(N, M, L);
    % utmp=pagetranspose(ncread(hisfile, 'u'));
    % vtmp=pagetranspose(ncread(hisfile, 'v'));
    wtmp=pagetranspose(ncread(hisfile, 'w'));
    % zetatmp=pagetranspose(ncread(hisfile, 'zeta'));
    % rhotmp=pagetranspose(ncread(hisfile, 'rho'));

    % construct 3d rho mat from get_rho.m
    for k=1:N
        vlevel = k;
        % tmp = get_rho(hisfile, gridfile, tindex, vlevel, coef);
        % size(tmp)
	[~, ~, ~, rhotmp(k, :, :)] = get_rho(hisfile, gridfile, tindex, vlevel, coef);
    end
    % utmp = shiftdim(utmp, 2);
    % vtmp = shiftdim(vtmp, 2);
    wtmp = shiftdim(wtmp, 2);
    % rhotmp = shiftdim(rhotmp, 2);

    size(wtmp)
    size(rhotmp)
    % size(zetatmp)
    % ub = ub + utmp;
    % vb = vb + vtmp;
    wb = wb + wtmp;
    % zetab = zetab + zetatmp;
    rhob = rhob + rhotmp;
end
% ub = ub./Nt;
% vb = vb./Nt;
wb = wb./Nt;
rhob = rhob./Nt;
% zetab = zetatmp./Nt;
%-----------------------------
%
% Save
%
% save data into a netcdf file
filename = strcat('nesm_2019_2020_5km_rho_w_annual_mean_3d_nt_',string(indxRange(1)), '_', ...
    string(indxRange(Nt)), '.nc');
ncid = netcdf.create(filename,'CLOBBER');
x_len = netcdf.defDim(ncid, 'Nx', Nx);
y_len = netcdf.defDim(ncid, 'Ny', Ny);
u_x_len = netcdf.defDim(ncid, 'Lu', Lu);
v_y_len = netcdf.defDim(ncid, 'Mv', Mv);
z_len = netcdf.defDim(ncid, 'Nz', N);
one = netcdf.defDim(ncid, 'one', 1);
netcdf.close(ncid);

% define variables and attributes
ncid = netcdf.open(filename,'WRITE');
netcdf.reDef(ncid);

% u_varid = netcdf.defVar(ncid, 'ub', 'double', [z_len y_len u_x_len]);
% netcdf.putAtt(ncid, u_varid, 'description', 'annual mean of u');
% netcdf.putAtt(ncid, u_varid, 'units', 'ms^-1');
% netcdf.putAtt(ncid, u_varid, 'array dimensions', size(ub));

% v_varid = netcdf.defVar(ncid, 'vb', 'double', [z_len v_y_len x_len]);
% netcdf.putAtt(ncid, v_varid, 'description', 'annual mean of v');
% netcdf.putAtt(ncid, v_varid, 'units', 'ms^-1');
% netcdf.putAtt(ncid, v_varid, 'array dimensions', size(vb));

w_varid = netcdf.defVar(ncid, 'wb', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, w_varid, 'description', 'annual mean of w');
netcdf.putAtt(ncid, w_varid, 'units', 'ms^-1');
netcdf.putAtt(ncid, w_varid, 'array dimensions', size(wb));

rho_varid = netcdf.defVar(ncid, 'rhob', 'double', [z_len y_len x_len]);
netcdf.putAtt(ncid, rho_varid, 'description', 'annual mean of rho');
netcdf.putAtt(ncid, rho_varid, 'units', 'kgm^-3');
netcdf.putAtt(ncid, rho_varid, 'array dimensions', size(rhob));

% zeta_varid = netcdf.defVar(ncid, 'zetab', 'double', [y_len x_len]);
% netcdf.putAtt(ncid, zeta_varid, 'description', 'annual mean of zeta');
% netcdf.putAtt(ncid, zeta_varid, 'units', 'm');
% netcdf.putAtt(ncid, zeta_varid, 'array dimensions', size(zetab));
% close define mode
netcdf.endDef(ncid);
%%
% put values
% netcdf.putVar(ncid, u_varid, ub); % start from 0
% netcdf.putVar(ncid, v_varid, vb); % start from 0
netcdf.putVar(ncid, w_varid, wb); % start from 0
% netcdf.putVar(ncid, zeta_varid, zetab); % start from 0
netcdf.putVar(ncid, rho_varid, rhob); % start from 0
% close netcdf file
netcdf.close(ncid);
disp('Done creating annual mean!')
