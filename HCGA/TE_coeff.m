function [coeff_TE] = TE_coeff(ksm,kam,a,s,kx0,n)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
u=kam.*a/2;
v=ksm.*s/2;
aTE=ksm.*sin(u)+kam.*cos(u).*sin(2*v)*exp(1i*kx0)+ksm.*sin(u).*cos(2*v)*exp(1i*kx0)
bTE=ksm.*cos(u)+kam.*sin(u).*sin(2*v).*exp(1i*kx0)-ksm.*cos(u).*cos(2*v).*exp(1i*kx0)
cTE=kam.*sin(v).*(exp(1i*kx0)+cos(2*u))+2*ksm.*cos(u).*sin(u).*cos(v)
dTE=-kam.*cos(v).*(exp(1i*kx0)-cos(2*u))-2*ksm.*cos(u).*sin(u).*sin(v)
coeff_TE=[aTE bTE cTE dTE];
% coeff_TE=coeff_TE./max(max(abs(coeff_TE)));
% coeff_TE
end

 