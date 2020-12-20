function [Hy] = plot_Region_TM(coeff,eigenmodes,kx0,l,ita,n,tg,lambda,difforder,meshdensity)

kss=eigenmodes;
eigenvalues=-sqrt(n^2/lambda^2-kss.^2/l^2);
m=size(eigenvalues,1);
for ind=1:m
    if(isreal(eigenvalues(ind)))
        eigenvalues(ind)=-eigenvalues(ind);
    end
end

beta=2*pi*eigenvalues;
a=(1-ita)*l;                            %thickness of air
s=(ita)*l;                              % thickness of bar
kam=sqrt(4*pi^2/lambda^2-beta.^2);      %x-wavevector in air region-II
ksm=sqrt(4*pi^2/lambda^2*n^2-beta.^2);  %x-wavevector in bar region-II
am=coeff(:,1);
bm=coeff(:,2);
cm=coeff(:,3);
dm=coeff(:,4);
x=linspace(0,3*l,round(3*l*meshdensity));
Hy1=(am.*cos(kam.*(mod(x,l)-a/2))+bm.*sin(kam.*(mod(x,l)-a/2))).*(exp(-2*i*kx0*l).*(((a+2*l)>x)&(x>2*l))+exp(-1*i*kx0*l).*(((a+l)>x)&(x>l))+((0<x)&(x<a)));
Hy2=(cm.*cos(ksm.*(mod(x,l)-(2*a+s)/2))+dm.*sin(ksm.*(mod(x,l)-(2*a+s)/2))).*((exp(-2*i*kx0*l).*(((3*l)>x)&(x>(a+2*l)))+exp(-1*i*kx0*l).*(((2*l)>x)&(x>(l+a)))+((x>a)&(x<l))));
Hy=Hy1+Hy2;
Ex=(beta/2/pi*lambda).*(Hy1+Hy2/n^2);
Hy(:,1)=[];
Hy(:,end)=[];
Ex(:,1)=[];
Ex(:,end)=[];
x(1)=[];
x(end)=[];
end