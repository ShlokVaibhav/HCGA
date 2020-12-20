function [coeff_TE] = TE_coeff(ksm,kam,a,s,kx0,n,l)
u=kam.*a/2;
v=ksm.*s/2;
aTE=ksm.*sin(u)+kam.*cos(u).*sin(2*v)*exp(1i*kx0*l)+ksm.*sin(u).*cos(2*v)*exp(1i*kx0*l);
bTE=ksm.*cos(u)+kam.*sin(u).*sin(2*v).*exp(1i*kx0*l)-ksm.*cos(u).*cos(2*v).*exp(1i*kx0*l);
cTE=kam.*sin(v).*(exp(1i*kx0*l)+cos(2*u))+2*ksm.*cos(u).*sin(u).*cos(v);
dTE=-kam.*cos(v).*(exp(1i*kx0*l)-cos(2*u))-2*ksm.*cos(u).*sin(u).*sin(v);
coeff_TE=[aTE bTE cTE dTE];
end

 