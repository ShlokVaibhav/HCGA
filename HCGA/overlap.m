function [Hnm,Enm] = overlap(kam,ksm,kx0,n1,cft,a,s,kxn,n,eigenvalues,gamman,type,l)
%% calculating the overlap integral H and E
m=size(kam,1);     %no. of modes
%% Now all the matrices are n*m, we can now vectorise the calculation for the overlap integral
Hnm=zeros(2*n1+1,m);
Hnm21=zeros(2*n1+1,m);       %n*m matrix
Hnm22=zeros(2*n1+1,m);       %n*m matrix
Enm=zeros(2*n1+1,m);         
for aa=1:2*n1+1              %for each n
    for bb=1:m               %for each (n,m)
        xn=kxn(aa);          %nth Diffracted wave's x-e=wavevector
        kam2=kam(bb);        %x-wavevector of mth mode in air
        ksm2=ksm(bb);        %x-wavevector of mth mode in bar
        Hnm21(aa,bb)=1./l./(xn^2-kam2^2).*(-(1+exp(1i*xn*a)).*sin(kam2*a/2).*(1i*cft(bb,2)*xn+cft(bb,1)*kam2)+(-1+exp(1i*xn*a)).*cos(kam2*a/2).*(-1i*cft(bb,1)*xn+cft(bb,2)*kam2));
        Hnm22(aa,bb)=1./l./(xn^2-ksm2^2)*(-(exp(1i*xn*a)+exp(1i*xn*l))*sin(ksm2*s/2)*(1i*cft(bb,4)*xn+cft(bb,3)*ksm2)-(exp(1i*xn*a)-exp(1i*xn*l))*cos(ksm2*s/2)*(cft(bb,4)*ksm2-1i*cft(bb,3)*xn));
        Hnm(aa,bb)=Hnm21(aa,bb)+Hnm22(aa,bb);
        if(type=='TM')
            Enm(aa,bb)=(Hnm21(aa,bb)+1/n^2*Hnm22(aa,bb))*eigenvalues(bb)/(gamman(aa))*2*pi;
        else
            Enm(aa,bb)=(Hnm21(aa,bb)+Hnm22(aa,bb))*(gamman(aa))/eigenvalues(bb)/2/pi;
        end
    end 
    
end
end


    

