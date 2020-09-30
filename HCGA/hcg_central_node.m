function [r_temp t_temp]=hcg_matrix2(gp,lambda,theta,ita, n_sub,n,tg,eigen,r_order)

% fignum=1;   % counts no. of figures
% timepass='N';
% meshdensity=10000;
% close all;
m1=size(eigen,1);
type='TM';                            
k=2*pi/lambda;                          %magnitude of incident wavevector
l=gp;                                    % The GP
theta=theta*pi/180;                     %angle converted to radians
difforder=round((m1-1)/2);                     
kx0=k*sin(theta);                       %x-component of incident wavevector
diffarray=-difforder:difforder;         %diffraction orders (difforder positive, difforder negative and reflection)
kxn=kx0+2*pi*diffarray/l;               %row of x-wavevector of all the diffraction orders
kxn=kxn.';                              %column of x-wavevectors

if(nargin==8)
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
tmain=HCG(gp,lambda,eigen(:,1),theta,ita(1),n(1),tg(1),type,m1);
for i=2:size(n,1)
 tmain=tmain*HCG(gp,lambda,eigen(:,i),theta,ita(i),n(i),tg(i),type,m1);
end

h_sub=eye(m1);
e_sub=eye(m1).*gamman_sub./gamman/(n_sub)^2;
T_I_sub=[(h_sub+e_sub)/2 (e_sub-h_sub)/2; (e_sub-h_sub)/2 (h_sub+e_sub)/2];
tmain=tmain*T_I_sub;

T11=tmain(1:2*difforder+1,1:2*difforder+1);
T12=tmain(1:1+2*difforder,2*difforder+2:4*difforder+2);
T21=tmain(2*difforder+2:4*difforder+2,1:1+2*difforder);
T22=tmain(2*difforder+2:4*difforder+2,2*difforder+2:4*difforder+2);
aout=inv(T11)*(ain-T12*bout);

bin=T21*aout;
T=abs(gamman_sub./(gamman(difforder+1)).*((abs(aout)).^2))./(n_sub)^2.*(gamman_sub==real(gamman_sub));
R=abs(gamman./(gamman(difforder+1)).*(abs(bin)).^2).*(gamman==real(gamman));
r_temp=sum(R);
t_temp=sum(T);

if(size(r_order,1)>0)
    r_order=r_order+difforder+1;
    for i=r_order
        r_temp=[r_temp;R(i)];
        t_temp=[t_temp;T(i)];
    end
end

end