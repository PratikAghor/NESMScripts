function Tadv=croco_vert_advection_aghor(maskr,omegaXomn_w,temp)
%
% compute Flux term in z direction
%
Vadv=croco_vert_trcflux(maskr,temp,omegaXomn_w);
%
Tadv=-Vadv;
%
return
