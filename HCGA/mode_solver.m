function [c]=get_modes(n,theta,ita,inn,num)
    c=[];
    x=inn;
    boundary=sqrt(inn.^2*(n^2-1)*ita^2);
    f=@(y)(((y*2*pi/ita).*(2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita)./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).*(1-ita)).*cos((y*2*pi/ita)*ita))+sinh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita)*(1-ita)).*sin((y*2*pi/ita)*ita).*(-(2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).^2+(y*2*pi/ita).^2./n^4))/(1+abs(cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)*ita^2)/ita).*(1-ita)))));
    g=@(y)(((y*2*pi/ita).*(2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita)./n^2.*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cos((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita).*(1-ita)).*cos((y*2*pi/ita)*ita))+sin((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita)*(1-ita)).*sin((y*2*pi/ita)*ita).*((2*pi*sqrt(y.^2-x.^2*(n^2-1)*ita^2)/ita).^2+(y*2*pi/ita).^2./n^4)));
    begin=0;
    if(boundary<0.1)
        c=[c FindRealRoots(f,0,boundary,500)];
    else
        while(begin<count)
        c=[c; FindRealRoots(f,begin,max(begin+1,count),100)];
        begin=begin+1;
        end
    end
    begin=count+0.0001;
    nummm=size(c,1);
    inn=2;
    iiii=1;
    for iii=1:nummm
        if(c(iiii)<0)
            c(iiii)=[];
            iiii=iiii-1;
        end
        iiii=iiii+1;
    end
    if (abs(c(1))<0.005)
        c(1)=[];
    end
            m=size(c,1);

    while inn<=m
        if(((c(inn)-c(inn-1)<0.01)||(abs(c(inn)))<0.01))
            c(inn)=[];
            m=m-1;
            inn=inn-1;
        end
        inn=inn+1;
    end
    
while(size(c,1)<num+2)
    c=[c;FindRealRoots(g,begin,begin+1,500)];
    
    m=size(c,1);
    begin=begin+1;
    inn=1;
    while inn<=m
    if(abs((c(inn)-count))<0.001)
       c(inn)=[];
       m=m-1;
       inn=inn-1;
   end
   inn=inn+1;
   end
    inn=2;
    while inn<=m
        if((c(inn)-c(inn-1))<0.001)
            c(inn)=[];
            m=m-1;
            inn=inn-1;
        end
        inn=inn+1;
    end


end

%      inn=2;
%         while inn<=m
%             if((c(inn)-c(inn-1))<0.01)
%                 c(inn)=[];
%                 m=m-1;
%                 inn=inn-1;
%             end
%             inn=inn+1;
%         end
%             

c=c(1:num);
end
