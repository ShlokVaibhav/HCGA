function [mid]=selfmade_fsolve_TM(n,theta,ita,lambda)

x=lambda;
val=[];
points=5e5;
countnn=sqrt(lambda^2*(n^2-1));
    ff=@(y)(real((1+sign(countnn-y)).*((y.*(sqrt(-y.^2+x.^2*(n^2-1)))./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sinh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)))*(1-ita)).*sin((y*2*pi)*ita).*(-(-y.^2+x.^2*(n^2-1))+y.^2./n^4))/(1+abs(cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita)))))));
    gg=@(y)((1+sign(y-countnn)).*(y.*(sqrt(y.^2-x.^2*(n^2-1)))./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cos((2*pi*sqrt(y.^2-x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sin((2*pi*sqrt(y.^2-x.^2*(n^2-1)))*(1-ita)).*sin(y*2*pi*ita).*((sqrt(y.^2-x.^2*(n^2-1))).^2+(y).^2./n^4)));
%     h=@(y)(real((1+sign(count-y)).*((y*2*pi/ita).*(2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita)./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).*(1-ita)).*cos((y*2*pi/ita)*ita))+sinh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita)*(1-ita)).*sin((y*2*pi/ita)*ita).*(-(2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).^2+(y*2*pi/ita).^2./n^4))/(1+abs(cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).*(1-ita))))+(1+sign(y-count)).*((y*2*pi/ita).*(2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita)./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cos((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita).*(1-ita)).*cos((y*2*pi/ita)*ita))+sin((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita)*(1-ita)).*sin((y*2*pi/ita)*ita).*((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita).^2+(y*2*pi/ita).^2./n^4))));
    h=@(y)(ff(y)+gg(y));

% s=linspace(countnn-0.005,countnn+0.005,points);
% ii=zeros(1000,1);
% for i=1:1000
%     ii(i)=h(countnn-0.005+0.01*(i-1)/1000);
% end
% plot(linspace(countnn-0.005,countnn+0.005,1000),ii,linspace(countnn-0.005,countnn+0.005,1000),0*ii)
low=countnn-0.0050;
high=countnn+0.0050;
mid=[];

if((h(low)*h(countnn-5e-5))<0)
    high=countnn-5e-5;
    mid=(low+high)/2;
    while(abs(low-high)>1e-6)
        if(h(mid)*h(high)<0)
            low=mid;
        elseif(h(mid)*h(low)<0)
            high=mid;
        else
            low=high;
        end
        mid=(low+high)/2;
    end
end
high=countnn+0.005;

    if((h(high)*h(countnn+5e-5))<0)

    low=countnn+5e-5;
    mid1=(low+high)/2;

    while(abs(low-high)>1e-6)
        if(h(mid1)*h(high)<0)
            low=mid1;
        elseif(h(mid1)*h(low)<0)
            high=mid1;
        else
            low=high;
        end
        mid1=(low+high)/2;
    end
    mid=[mid;mid1];
    end
if(h(high)*h(low)<0)
    if(size(mid,1)==0)
        mid=countnn;
    end
end

end