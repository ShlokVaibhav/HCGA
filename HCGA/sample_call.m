%% This is a sample script to show how to call the functions, the RCWA
%  function can be called to verify the HCGA results as well
%  RCWA can be called only in TM mode
close all
tic;
grating_period=1;
theta=0;
ita=0.4;
width1=0.2;
theta=10;
num=25;

n_substrate=1;
n1=3.48;
%%
r_a_1=zeros(1000,1);
r_a_2=r_a_1;
t_a_2=r_a_1;
t_a_1=zeros(1000,1);
r_r_1=r_a_1;
t_r_1=r_a_1;
if(~mod(num,2))
    num=num+1;
end
%% To find and plot reflectivity and compare with RCWA result
for i=1:1000
    wavelength=1+i*0.002;
[r_a_1(i), t_a_1(i)]=hcg_central_node(grating_period,wavelength, theta, [ita], n_substrate,[n1],[width1],num,'TM','N');
[r_r_1(i), t_r_1(i)]=RTA_1d_tm(grating_period,1,1,(n_substrate)^2+0i,...
                       [n1^2+0i],[n1^2+0i],[1+0i],...
                       [1+0i],[1-ita],[1;width1;1],35,2*pi/wavelength,2*pi/wavelength*sin(theta*pi/180));

end
plot(linspace(1.002,3,1000), r_a_1,linspace(1.002,3,1000), r_r_1,'--', 'LineWidth', 2)
set(0,'DefaultTextInterpreter','Latex')
xlabel '$\lambda$ ($\mu$ m)'
ylabel 'Reflectivity'
legend 'HCGA reflectivity' 'RCWA Reflectivity'

%% To find and plot field profile
wavelength=2;
[r_a_1(i), t_a_1(i)]=hcg_central_node(grating_period,wavelength, theta, [ita], n_substrate,[n1; 2],[width1; width1],num,'TE','Y');


toc;
