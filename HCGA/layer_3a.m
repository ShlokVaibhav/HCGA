
% cd('D:\\Documents\\Matlab\\HCG_trial')

%% Fill details in this section
%This is for 3 layered structure
warning('off');
grating_size=1;   
width_size=1;
duty_size=1;
grating=linspace(2,2,grating_size);
width_array=linspace(2.5,2.5,width_size)
duty_array=linspace(0.7694,0.7694,duty_size)
wavelength_size=100;
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
try
 
r_r_off=zeros(grating_size,width_size,duty_size,wavelength_size);
t_r_off=r_r_off;
r_r_on=r_r_off;
t_r_on=r_r_off;
r_a_off=zeros(grating_size,width_size,duty_size,wavelength_size);
r_a_on=r_a_off;
t_a_off=r_a_off;
t_a_on=r_a_off;
em1=[];
em2=[];
em3=[];
step=0;
steps=grating_size*width_size*wavelength_size*duty_size;
fff = waitbar(0,'1','Name','Reflectivity calculating...',...
               'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
%      

for grating_iter=1:grating_size
    tic
    grating_period=grating(grating_iter);
    delta=delta/10;
%     num=2*grating_iter-1;
    for width_iter=1:width_size
        wide=width_array(width_iter);
        for ita_iter=1:duty_size
            ita=duty_array(ita_iter);
            if(1==1)
             if(1==1)
                for wavelength_iter=1:wavelength_size
%                      lower_wavelength=peak_wavelength_2900nm(width_iter,ita_iter)-50e-6;
%                      upper_wavelength=lower_wavelength+200e-6;

                    wavelength_array=linspace(lower_wavelength,upper_wavelength,wavelength_size);
                    step=step+1;
                    wavelength = wavelength_array(wavelength_iter);
                    width_GaN_bottom = wide;
                    waitbar(step/steps,fff, sprintf('Processing %d of %d...',step,steps))
                    if getappdata(fff,'canceling')
                        break
                    end

                    width_AlGaN=0.025;
                    width_2DEG=0.01;
                    width_GaN_mid=0.1;
                  
                    theta=0;

                   if(~mod(num,2))
                        num=num+1;
                   end
                   
                   if(compute_analytical=='Y')
                       em1=[em1 get_modes(n_AlGaN,theta,ita,grating_period/wavelength,num)];
                       em2=[em2 get_modes(n_2DEG_on,theta,ita,grating_period/wavelength,num)];
                       em3=[em3 get_modes(n_GaN,theta,ita,grating_period/wavelength,num)];
                       [r_a_on(grating_iter,width_iter,ita_iter,wavelength_iter,:) t_a_on(grating_iter,width_iter,ita_iter,wavelength_iter,:)]=hcg_central_node(grating_period,wavelength, theta, [ita;ita;ita], n_substrate,[n_AlGaN;n_2DEG_on;n_GaN],[width_AlGaN;width_2DEG;width_GaN_bottom],[em1(1:end,end) em2(1:end,end) em3(1:end,end)]);
%                        [r_a_off(grating_iter,width_iter,ita_iter,wavelength_iter,:) t_a_off(grating_iter,width_iter,ita_iter,wavelength_iter,:)]=hcg_central_node(grating_period,wavelength, theta, [ita;ita;ita], n_substrate,[n_AlGaN;n_2DEG_off;n_GaN],[width_AlGaN;width_2DEG;width_GaN_bottom],[em1(1:end,end) em3(1:end,end) em3(1:end,end)]);

                   end
                   if(compute_RCWA=='Y')
                       [r_r_on(grating_iter,width_iter,ita_iter,wavelength_iter) t_r_on(grating_iter,width_iter,ita_iter,wavelength_iter)]=RTA_1d_tm(grating_period,3,1,(n_substrate)^2+0i,...
                       [n_AlGaN^2+0i;n_2DEG_on^2+0i;n_GaN^2+0i],[n_AlGaN^2+0i;n_2DEG_on^2+0i;n_GaN^2+0i],[1+0i;1+0i;1+0i],...
                       [1+0i;1+0i;1+0i],[1-ita;1-ita;1-ita],[1;  width_AlGaN;width_2DEG;width_GaN_bottom;1],300,2*pi/wavelength,2*pi/wavelength*sin(theta*pi/180));
%                       [r_r_off(grating_iter,width_iter,ita_iter,wavelength_iter) t_r_off(grating_iter,width_iter,ita_iter,wavelength_iter)]=RTA_1d_tm(grating_period,3,1,(n_substrate)^2+0i,...
%                        [n_AlGaN^2+0i;n_2DEG_off^2+0i;n_GaN^2+0i],[n_AlGaN^2+0i;n_2DEG_off^2+0i;n_GaN^2+0i],[1+0i;1+0i;1+0i],...
%                        [1+0i;1+0i;1+0i],[1-ita;1-ita;1-ita],[1;  width_AlGaN;width_2DEG;width_GaN_bottom;1],300,2*pi/wavelength,2*pi/wavelength*sin(theta*pi/180));    
                   end
                end
            else
                    em1=[em1 zeros(num,wavelength_size)];
                    em2=[em2 zeros(num,wavelength_size)];
                    em3=[em3 zeros(num,wavelength_size)];

                    step=step+wavelength_size;
                    waitbar(step/steps,fff, sprintf('Processing %d of %d... but skipping',step,steps))
                    if getappdata(fff,'canceling')
                        break
                    end
               
            end

        end 
%     time_analytical(grating_iter,width_iter)=toc;
    step
    end
    end
end
toc
 delete(fff)
beep
% 
% if(compute_RCWA=='Y')
%         path=pwd;
%     cd 'D:\Documents\MATLAB\HCG_trial\Analytical 2.9um'
%   save(['r_r_on_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'r_r_on')
%   save(['r_r_off_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'r_r_off')
%   cd(path);
% 
% end
% if(compute_analytical=='Y')
%     path=pwd;
%   cd 'D:\Documents\MATLAB\HCG_trial\Analytical 2.9um'
%   save(['r_a_on_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'r_a_on')
%   save(['r_a_off_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'r_a_off')
%   save(['t_a_on_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'t_a_on')
%    save(['t_a_off_width_',num2str(width_array(1)),'um_',num2str(width_array(end)),'_um_',num2str(width_size),'_samples_grating_',num2str(grating),'_um_duty_',num2str(duty_array(1)),'_',num2str(duty_array(end)),'_',num2str(duty_size),'_samples_',num2str(wavelength_size),'_points_',num2str(lower_wavelength),'_',num2str(upper_wavelength),'_nm_.mat'] ,'t_a_off')
% 
%   cd(path);
% end



catch ME
    delete(fff);
    beep
    rethrow(ME);
end
