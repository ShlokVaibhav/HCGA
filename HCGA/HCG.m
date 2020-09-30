function [T_region, coefficients] = HCG(gp,lambda,eigenmodes,theta,ita,n,tg,type,m)
l=gp;
k=2*pi/lambda;
difforder=round((m-1)/2);                %no. of positive diffracted waves to be considered
kx0=k*sin(theta);                        %x-component of incident wavevector
diffarray=-difforder:difforder;          %diffraction orders (positive, negative and reflection)
kxn=kx0+2*pi*diffarray/l;                %row ofx-wavevector of diffraction orders
kxn=kxn.';                               %column of x-wavevectors
gamman=-sqrt(k^2-kxn.^2);                %z-component of wavevectors in air
for nn=1:2*difforder+1
    if(isreal((gamman(nn))))
        gamman(nn)=sqrt(k^2-kxn(nn)^2);  %real z-wavevector of diffracted rays
    end
end


if(ita==1)                               % homogenouos layer
    beta=-sqrt(n^2*k^2-kxn.^2);          % z-component of wavevector in the substrate
    for nn=1:2*difforder+1               
        if(isreal((beta(nn))))
            beta(nn)=sqrt(n^2*k^2-kxn(nn)^2);  %z-wavevector of real diffracted rays be positive
        end
    end
    hoverlap=eye(2*difforder+1);
    eoverlap=eye(2*difforder+1).*(gamman./beta)*n^2;
    eoverlap=pinv(eoverlap);
    T_I_HCG=[(hoverlap+eoverlap)/2 (-hoverlap+eoverlap)/2 ;(-hoverlap+eoverlap)/2 (hoverlap+eoverlap)/2];
    T_HCG_III=pinv(T_I_HCG);
else
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
    kam=sqrt(k^2-beta.^2);   %x-wavevector in air region-II
    ksm=sqrt(k^2*n^2-beta.^2);    %x-wavevector in bar region-II
    if(type=='TE')
        coefficients=TE_coeff(ksm,kam,a,s,kx0,n);     %for TE mode
    elseif(type=='TM')
        coefficients=TM_coeff(ksm,kam,a,s,kx0,n,l);     %for TM mode
    end
%   H2=plot_Region2(coefficients(:,1:4),kam,ksm,a,s,kx0,l,40,100,beta,tg);
    [hoverlap, eoverlap]=overlap2(kam,ksm,kx0,difforder,coefficients(:,1:4),a,s,kxn,n,eigenvalues,gamman,type,l);
    T_I_HCG=[(hoverlap+eoverlap)/2 (-hoverlap+eoverlap)/2 ;(-hoverlap+eoverlap)/2 (hoverlap+eoverlap)/2];
%     [hoverlap2, eoverlap2]=overlap23(kam,ksm,kx0,difforder,coefficients(:,1:5),a,s,kxn,n,eigenvalues,gamman,type,l);
%     T_HCG_III=[(hoverlap2+eoverlap2)/2 (-hoverlap2+eoverlap2)/2 ;(-hoverlap2+eoverlap2)/2 (hoverlap2+eoverlap2)/2];
    T_HCG_III=inv(T_I_HCG);
    
end
m=size(beta,1);
psi=diag(exp(-1i*tg*beta));
T_HCG=[inv(psi) zeros(m,m); zeros(m,m) psi];
cutoff=1e3;
for i=1:2*m
        if(abs(T_HCG(i,i))<1/cutoff)
            T_HCG(i,i)=1/cutoff;
        elseif(abs(T_HCG(i,i))>cutoff)
            T_HCG(i,i)=cutoff;
        end
end
                                            
T_region=T_I_HCG*T_HCG*T_HCG_III;
end

