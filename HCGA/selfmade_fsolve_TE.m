function [mid]=selfmade_fsolve_TE(n,theta,ita,lambda_inv)
    x=lambda_inv;
    val=[];
    points=5e5;                         % Fineness of grid
    countnn=sqrt(lambda_inv^2*(n^2-1)); % Cutoff (where k_a=0)
    ff=@(y)(real((1+sign(countnn-y)).*((y.*(sqrt(-y.^2+x.^2*(n^2-1))).*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sinh((2*pi*sqrt(-y.^2+x.^2*(n^2-1)))*(1-ita)).*sin((y*2*pi)*ita).*(-(-y.^2+x.^2*(n^2-1))+y.^2))/(1+abs(cosh((2*pi*sqrt(-y.^2+x.^2*(n^2-1))).*(1-ita)))))));
    %k_a imaginary
    gg=@(y)((1+sign(y-countnn)).*(y.*(sqrt(y.^2-x.^2*(n^2-1))).*(2*cos(2*pi.*x.*sin(theta*pi/180))-2*cos((2*pi*sqrt(y.^2-x.^2*(n^2-1))).*(1-ita)).*cos(y*2*pi*ita))+sin((2*pi*sqrt(y.^2-x.^2*(n^2-1)))*(1-ita)).*sin(y*2*pi*ita).*((sqrt(y.^2-x.^2*(n^2-1))).^2+(y).^2)));
    %k_a real
    h=@(y)(ff(y)+gg(y));
    %linear combination valid in both cases
    low=countnn-0.0050;                 % lower search boundary 
    high=countnn+0.0050;                % high search boundary
    mid=[];                             % push root candidates into it

    if((h(low)*h(countnn-5e-5))<0)      %If a root is present in left half
        high=countnn-5e-5;              % Do standard binary search
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

    if((h(high)*h(countnn+5e-5))<0)     % If root is present in right-half
    low=countnn+5e-5;                   % Do standard binary search in right-half
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
        if(size(mid,1)==0)            %If no root found in binary search but there is
            mid=countnn;              % one, then it is very close to cutoff
        end                           % in this case, assign cutoff as root
    end

end