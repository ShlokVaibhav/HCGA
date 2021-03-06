function [r_temp, t_temp]=hcg_central_node(gp,lambda,theta,ita, n_sub,n,tg,num_eigen,type,plot_enable,r_order)
% gp - Grating period
% lambda - wavelength
% theta - in degress
% ita - duty cycle (of bar)
% 
% ALL REFRATIVE INDICES MUST BE REAL
% n_sub - refractive index of substrate
% n - m X 1 vector , bar refractive index for each layer
% tg - m X 1 vector , thickness of each layer
% num_eigen - no. of modes used, must be even, 25 is good
% type  - 'TE' or 'TM'
% plot_enable - 'Y' to plot field profile, 'N' to not
% r_order - provide if reflectivity for a particular mode is needed, beta stage

mesh_density=1000;
eigen=[];
for i=1:size(tg,1)
    if(type=='TE')
        eigen=[eigen get_modes_TE(n(i),theta,ita,gp/lambda,num_eigen)];
    else
        eigen=[eigen get_modes_TM(n(i),theta,ita,gp/lambda,num_eigen)];
    end
end
cache = {};
m1=num_eigen;
k=2*pi/lambda;                          %magnitude of incident wavevector
l=gp;                                   % The GP
theta=theta*pi/180;                     %angle converted to radians
difforder=round((m1-1)/2);                     
kx0=k*sin(theta);                       %x-component of incident wavevector
diffarray=-difforder:difforder;         %diffraction orders (difforder positive, difforder negative and reflection)
kxn=kx0+2*pi*diffarray/l;               %row of x-wavevector of all the diffraction orders
kxn=kxn.';                              %column of x-wavevectors

if(nargin==10)
    r_order=[];
end
gamman=-sqrt(k^2-kxn.^2);               %z-component of wavevectors in Air containing Incident wave
gamman_sub=-sqrt((n_sub)^2*k^2-kxn.^2); %z-component of wavevectors in Substrate
for nn=1:2*difforder+1
    if(isreal((gamman(nn))))
        gamman(nn)=sqrt(k^2-kxn(nn)^2);          %real z-wavevector in air should have positive sign
    end

    if(isreal((gamman_sub(nn))))
        gamman_sub(nn)=sqrt(n_sub^2*k^2-kxn(nn)^2); %real z-wavevector in substrate should have positive sign
    end
end
ain2=zeros(difforder,1);
ain=[ain2 ;1;ain2];

bout=zeros(2*difforder+1,1);
tmain=eye(4*difforder+2);
for i=1:size(n,1)
    [a,b,c,d,e]=HCG(gp,lambda,eigen(:,i),theta,ita,n(i),tg(i),type,m1,difforder);
    cache{end+1}=a;
    cache{end+1}=b;
    cache{end+1}=c;
    cache{end+1}=d;
    cache{end+1}=e;
    tmain=tmain*a;
end

h_sub=eye(2*difforder+1);
if(type=='TE')
    e_sub=eye(2*difforder+1).*gamman./gamman_sub;
else
    e_sub=eye(m1).*gamman_sub./gamman/(n_sub)^2;
end
T_I_sub=[(h_sub+e_sub)/2 (e_sub-h_sub)/2; (e_sub-h_sub)/2 (h_sub+e_sub)/2];
tmain=tmain*T_I_sub;
T11=tmain(1:2*difforder+1,1:2*difforder+1);
T12=tmain(1:1+2*difforder,2*difforder+2:4*difforder+2);
T21=tmain(2*difforder+2:4*difforder+2,1:1+2*difforder);
T22=tmain(2*difforder+2:4*difforder+2,2*difforder+2:4*difforder+2);

aout=inv(T11)*(ain-T12*bout);
bin=T21*aout;
cft=[aout; bout];
length=0;
pointss=0;
field_HCG=[];
if(plot_enable=='Y')
for kk=1:size(tg,1)
     i=size(tg,1)+1-kk;
     matr = cache{4+5*(i-1)};
     beta = cache{3+5*(i-1)};
     count= cache{5+5*(i-1)};
     weights = cache{2+5*(i-1)};
     T_region=cache{1+5*(i-1)};
     cftnew=matr*cft;
     cft=T_region*cft;
     cft-[ain;bin];
     a=cftnew(1:m1);
     b=cftnew(m1+1:end);
     a((m1-count+1):end)=[];
     b((m1-count+1):end)=[];
     z=linspace(sum(tg)-length-tg(i), sum(tg)-length,round(tg(i)*mesh_density));
     pointss=pointss+size(z,2);
     beta((m1-count+1):end)=[];
     if(type=='TE')
         field_z=(a).*exp(-1j*beta*(z-sum(tg)+length))+(b).*exp(1j*beta*(z-sum(tg)+length));
     else
         field_z=(a).*exp(-1j*beta*(z-sum(tg)+length))-(b).*exp(1j*beta*(z-sum(tg)+length));
     end
     length=length+tg(i);
     if(type=='TE')
        field_x=plot_Region_TE(weights(1:(m1-count),:),eigen(1:(m1-count),i),kx0,l,ita,n(i),tg(i),lambda,difforder,mesh_density);
     else
        field_x=plot_Region_TM(weights(1:(m1-count),:),eigen(1:(m1-count),i),kx0,l,ita,n(i),tg(i),lambda,difforder,mesh_density);
     end
     field = zeros(size(field_z,2),size(field_x,2));
     for i=1:(m1-count)
        field=field+(field_z(i,:)).'*(field_x(i,:));
     end
        field_HCG=[field;field_HCG];
end
     ln=1;
     z1=linspace(-ln,0,round(ln*mesh_density));
     pointss=pointss+size(z1,2);
     if(type=='TE')
        field_z1=-k./gamman.*ain.*exp(-1j*gamman*z1)-k./gamman.*bin.*exp(1j*gamman*z1);
     else
        field_z1=ain.*exp(-1j*gamman*z1)-bin.*exp(1j*gamman*z1);
     end
     field_x2=linspace(0,3*l,round(3*l*mesh_density));
     field_x2(1)=[];
     field_x2(end)=[];
     field_x2=exp(-1j*kxn*field_x2);
     fieldn=zeros(size(field_z1,2),size(field_x2,2));
     for i=1:m1
          fieldn=fieldn+(field_z1(i,:)).'*(field_x2(i,:));
     end
     fieldn=[fieldn;field_HCG];
     ln1=0.3;
     z1=linspace(sum(tg),sum(tg)+0.2,round(ln1*mesh_density));
     pointss=pointss+size(z1,2);
     if(type=='TE')
     field_z1=-k./gamman.*aout.*exp(-1j*gamman*(z1-sum(tg)))-k./gamman.*bout.*exp(1j*gamman*(z1-sum(tg)));
     else
     field_z1=aout.*exp(-1j*gamman*(z1-sum(tg)))-bout.*exp(1j*gamman*(z1-sum(tg)));
     end
     field_x2=linspace(0,3*l,round(3*l*mesh_density));
     field_x2(1)=[];
     field_x2(end)=[];
     field_x2=exp(-1j*kxn*field_x2);
     field_x2=field_x2;
     fieldnn=zeros(size(field_z1,2),size(field_x2,2));
     for i=1:m1
          fieldnn=fieldnn+(field_z1(i,:)).'*(field_x2(i,:));
     end
     fieldn=[fieldn;fieldnn];
      figure
     xax=linspace(0,3*l,round(3*l*mesh_density));
     xax(1)=[];
     xax(end)=[];
     imagesc(xax,linspace(-ln,ln1+sum(tg), pointss),abs(fieldn).^2)
     colormap(jet)
     length=0;
     length_prev=0;
     for i=1:size(tg,1)
     length_prev=length_prev+length;
     length=tg(i)
     rectangle('Position',[(1-ita)*l,length_prev,ita*l,length],'EdgeColor','w','LineStyle','--','LineWidth',1)
     rectangle('Position',[(1-ita)*l+l,length_prev,ita*l,length],'EdgeColor','w','LineStyle','--','LineWidth',1)
     rectangle('Position',[(1-ita)*l+2*l,length_prev,ita*l,length],'EdgeColor','w','LineStyle','--','LineWidth',1)
     end
     set(0,'DefaultTextInterpreter','Latex')
     xlabel 'x ($\mu$m)'
     ylabel 'z ($\mu$m)'
     dim = [0.20 0.8 0.1 0.1];
     if(type=='TE')
        str = 'Analytical |E_y|^2';
     else
        str = 'Analytical |H_y|^2';
     end
     annotation('textbox',dim,'String',str,'FontSize',16)
     set(gca,'FontSize',24)
     colorbar
end
%% Reflectivity and Transmittivity module
if(type=='TE')
T=abs((gamman(difforder+1)./gamman_sub).*((abs(aout)).^2)).*(gamman_sub==real(gamman_sub));
R=abs((gamman(difforder+1))./gamman.*(abs(bin)).^2).*(gamman==real(gamman));
else
T=abs(gamman_sub./(gamman(difforder+1)).*((abs(aout)).^2))./(n_sub)^2.*(gamman_sub==real(gamman_sub));
R=abs(gamman./(gamman(difforder+1)).*(abs(bin)).^2).*(gamman==real(gamman));
end
r_temp=sum(R);
t_temp=sum(T);

if(size(r_order,1)>0)
    r_order=r_order+difforder+1;
    for i=r_order
        r_temp=[R(i)];
        t_temp=[T(i)];
    end
end
     
end