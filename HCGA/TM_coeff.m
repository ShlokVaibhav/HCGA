function [coeff_TM] = TM_coeff(ksm,kam,a,s,kx0,n,l)
u=kam.*a/2;
v=ksm.*s/2;
aTM=1./n^2.*ksm.*sin(u)+kam.*cos(u).*sin(2*v)*exp(1i*kx0*l)+1./n^2.*ksm.*sin(u).*cos(2*v)*exp(1i*kx0*l);
bTM=1./n^2.*ksm.*cos(u)+kam.*sin(u).*sin(2*v).*exp(1i*kx0*l)-1./n^2.*ksm.*cos(u).*cos(2*v).*exp(1i*kx0*l);
cTM=kam.*sin(v).*(exp(1i*kx0*l)+cos(2*u))+2./n^2.*ksm.*cos(u).*sin(u).*cos(v);
dTM=-kam.*cos(v).*(exp(1i*kx0*l)-cos(2*u))-2./n^2.*ksm.*cos(u).*sin(u).*sin(v);
coeff_TM=[aTM bTM cTM dTM];
end


 