function [c]=get_modes_TE(n,theta,ita,inn,num)
    x=inn;
    count=sqrt(inn.^2*(n^2-1));
    % cutoff where k_a = 0
    f=@(y)((1+sign(count-y))*((y.*(sqrt(-y.^2+x.^2*(n^2-1))).*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sinh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)))*(1-ita)).*sin((y*2*pi)*ita).*(-(sqrt(-y.^2+x.^2*(n^2-1))).^2+y.^2))/(1+abs(cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita))))));
    % when k_a is imaginary
    g=@(y)((1+sign(y-count))*(y.*(sqrt(y.^2-x.^2*(n^2-1))).*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cos((2*pi*sqrt(y.^2-x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sin((2*pi*sqrt(y.^2-x.^2*(n^2-1)))*(1-ita)).*sin(y*2*pi*ita).*((sqrt(y.^2-x.^2*(n^2-1))).^2+(y).^2)));
    % when k_a is real
    h=@(y)f(y)+g(y);
    
    begin=0;
    up=1;
    c=[];
    begin=0;
    if(count<0.1)
        c=[c FindRealRoots(f,0,count-0.01,500)];
    else
        while(begin<count-0.005)
            if(begin+1>count-0.005)
                c=[c; FindRealRoots(f,begin,min(begin+1,count-0.005),50)];
            else
                c=[c; FindRealRoots(f,begin,min(begin+1,count-0.005),50)];
            end
                begin=begin+1;
        end
    end
%      c=[c; FindRealRoots(f,count-0.001,count-0.000001,50)];

    
    begin=count+0.005;
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
    if(size(c,2)>0)
        if (abs(c(1))<0.01)
            c(1)=[];
        end
    end
    d=selfmade_fsolve_TE(n,theta,ita,x);
    c=[c;d];
    m=size(c,1);
     if(m>1)
         inn=2;
         eliminate=zeros(size(c));
         while inn<=m
             if((h((c(inn)+1e-4))*h(c(inn-1)-1e-4))<0)
                 eliminate(inn)=1;
                 c(inn);
                 h((c(inn)+1e-4));
                 h(c(inn-1)-1e-4);
             end
             inn=inn+1;
         end
         c(find(eliminate==1))=[];
     end
     m=size(c,1);
     iter=1;
     begin=count+0.005;
while(size(c,1)<num+2)
    if(begin==count+0.005)
        cc=FindRealRoots(g,begin,begin+0.5,50);
        if(size(cc,1)==0)
            cc=FindRealRoots(g,begin,begin+0.5,100);
        end
        c=[c;cc];
    else
        cc=FindRealRoots(g,begin,begin+0.5,50);
        if(size(cc,1)==0)
            cc=FindRealRoots(g,begin,begin+0.5,100);
        end
        c=[c;cc];
    end
    mm=size(c,1);
    begin=begin+0.5;
          inn=2;
         
         eliminate=zeros(size(c));
         while inn<=mm
             if((h((c(inn)-1e-4))*h(c(inn-1)+1e-4))<0)
                 eliminate(inn)=1;
                 c(inn-1);
                  c(inn);
                  (h((c(inn)-1e-4)));
                   h(c(inn-1)+1e-4);
%                    
             end
             inn=inn+1;
         end
         c(find(eliminate==1))=[];
     iter=iter+1;

end
  c=c(1:num);

end