function [coeff_TM] = TM_coeff(ksm,kam,a,s,kx0,n,l)
%%
%   Detailed explanation goes here
u=kam.*a/2;
v=ksm.*s/2;
aTM=1./n^2.*ksm.*sin(u)+kam.*cos(u).*sin(2*v)*exp(1i*kx0*l)+1./n^2.*ksm.*sin(u).*cos(2*v)*exp(1i*kx0*l);
bTM=1./n^2.*ksm.*cos(u)+kam.*sin(u).*sin(2*v).*exp(1i*kx0*l)-1./n^2.*ksm.*cos(u).*cos(2*v).*exp(1i*kx0*l);
cTM=kam.*sin(v).*(exp(1i*kx0*l)+cos(2*u))+2./n^2.*ksm.*cos(u).*sin(u).*cos(v);

dTM=-kam.*cos(v).*(exp(1i*kx0*l)-cos(2*u))-2./n^2.*ksm.*cos(u).*sin(u).*sin(v);
coeff_TM=[aTM bTM cTM dTM];

 size(max(coeff_TM,[],2));
aTM=coeff_TM(:,1);
bTM=coeff_TM(:,2);
cTM=coeff_TM(:,3);
dTM=coeff_TM(:,4);
c=zeros(size(kam));
for i=1:size(kam,1)
    if(isreal(kam(i)))
        c(i)=((abs(aTM(i))).^2+(abs(bTM(i))).^2).*a./2+((abs(cTM(i))).^2+(abs(dTM(i))).^2).*s./2./n^2+((abs(aTM(i))).^2-(abs(bTM(i))).^2).*sin(kam(i).*a)/2./kam(i)+((abs(cTM(i))).^2-(abs(dTM(i))).^2).*sin(ksm(i).*s)./2./ksm(i)./n^2;
    else
        c(i)=((abs(aTM(i))).^2+(abs(bTM(i))).^2).*sinh(abs(kam(i)).*a)./2./abs(kam(i))+((abs(cTM(i))).^2+(abs(dTM(i))).^2).*s/2./n^2+((abs(aTM(i))).^2-(abs(bTM(i))).^2).*a./2+((abs(cTM(i))).^2-(abs(dTM(i))).^2).*sin(ksm(i).*s)./2./ksm(i)./n^2;
    end
end
% 
% h=ones(size(kam,1),size(kam,1));
coeff_TM=[coeff_TM c];
% for i=1:size(kam,1)
%     for j=1:size(kam,1)
%         if(i~=j)
%             am1=aTM(i);
%             am2=conj(bTM(j));
%             bm1=aTM(i);
%             bm2=conj(bTM(j));
%             cm1=cTM(i);
%             cm2=conj(cTM(j));
%             dm1=dTM(i);
%             dm2=conj(dTM(j));
%             ka1=kam(i);
%             ka2=conj(kam(j));
%             ks1=conj(ksm(i));
%             ks2=conj(ksm(j));
%             h(i,j)=(am1*am2+bm1*bm2)*sin((ka1-ka2)*a/2)/(ka1-ka2)+(am1*am2-bm1*bm2)*sin((ka1+ka2)*a/2)/(ka1+ka2)+(cm1*cm2-dm1*dm2)*sin((ks1+ks2)*s/2)/(ks1+ks2)/n^2+(cm1*cm2+dm1*dm2)*sin((ks1-ks2)*s/2)/(ks1-ks2)/n^2;
%         else
%             h(i,j)=c(i);
%         end
%     end
end


 