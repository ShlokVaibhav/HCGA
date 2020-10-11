
% cd('D:\\Documents\\Matlab\\HCG_trial')

%% Fill details in this section
%This is for 3 layered structure
warning('off');
lower_wavelength=2;
upper_wavelength=4;
n_substrate=1;
n_AlGaN=2.15;
n_2DEG_on=2.3065;
n_2DEG_off=2.31;
n_GaN=2.31;
num=150;
compute_analytical='Y';
compute_RCWA='N';
delta=0;
%% 

% figure(1)
% rectangle('Position',[0,0,grating_period(1)*ita(1),width_GaN(1)*grating_period(1)],'FaceColor',[1 0 1],'EdgeColor','w')
% hold on
% rectangle('Position',[0,0+width_GaN*grating_period,grating_period*ita,width_2DEG*grating_period],'FaceColor',[0 1 1],'EdgeColor','w')
% rectangle('Position',[0,0+width_GaN*grating_period+width_2DEG,grating_period*ita,width_AlGaN*grating_period],'FaceColor',[1 1 0],'EdgeColor','w')
% rectangle('Position',[0+grating_period,0,grating_period*ita,width_GaN*grating_period],'FaceColor',[1 0 1],'EdgeColor','w')
% rectangle('Position',[0+grating_period,0+width_GaN*grating_period,grating_period*ita,width_2DEG*grating_period],'FaceColor',[0 1 1],'EdgeColor','w')
% rectangle('Position',[0+grating_period,0+width_GaN*grating_period+width_2DEG*grating_period,grating_period*ita,width_AlGaN*grating_period],'FaceColor',[1 1 0],'EdgeColor','w')
% axis([0 2*grating_period+0.1 0 (width_GaN+width_AlGaN+width_2DEG)*grating_period+0.1] )

                    wavelength = wavelength_array(wavelength_iter);
                    width_GaN_bottom = wide;
                    width_AlGaN=0.025;
                    width_2DEG=0.01;
                    width_GaN_mid=0.1;
                    theta=0;
                   if(~mod(num,2))
                        num=num+1;
                   end
                       em1=[em1 get_modes(n_AlGaN,theta,ita,grating_period/wavelength,num)];
                       em2=[em2 get_modes(n_2DEG_on,theta,ita,grating_period/wavelength,num)];
                       em3=[em3 get_modes(n_GaN,theta,ita,grating_period/wavelength,num)];
                       [r_a_on(grating_iter,width_iter,ita_iter,wavelength_iter,:) t_a_on(grating_iter,width_iter,ita_iter,wavelength_iter,:)]=hcg_central_node(grating_period,wavelength, theta, [ita;ita;ita], n_substrate,[n_AlGaN;n_2DEG_on;n_GaN],[width_AlGaN;width_2DEG;width_GaN_bottom],[em1(1:end,end) em2(1:end,end) em3(1:end,end)]);

                   
 