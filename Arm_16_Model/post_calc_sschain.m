% see if we can find a way to post-calculate the sschain from the chain
% results 
clc
clear
close all 

final_time = 0.5; 
time_interval2 = final_time/100;
time_int = (0:time_interval2:final_time)';

% Input_raw = dlmread('arm16_pert4_states.sto','\t',7,0);
input_struct = load('arm16_noise_kinematics.mat');
input_controls_ref = readmatrix('arm16_Tracking_p50_degroote_w75_cubed_v6controls.sto','FileType','text');
input_data = input_struct.new_noise_data; 

position = interp1(input_data(:,1),input_data(:,2),time_int);
velocity = interp1(input_data(:,1),input_data(:,3),time_int);
data.xdata = time_int;
data.ydata = [time_int,position,velocity];

load chain_results_20220115T023505

mu = [];
sig = [];

n_pools = size(chain,3);

queue = 1:1000:150000;

ss = zeros(length(queue),size(chain,3));
prior = zeros(length(queue),size(chain,3)); 
sschain = zeros(length(queue),size(chain,3)); 


tic 

for j = 1:n_pools
    for i = 1:length(queue)
        ss(i,j) = Arm16_SS(chain(queue(i),:,j),data);
        prior(i,j) = prior_act(chain(queue(i),:,j),mu,sig);    
        sschain(i,j) = ss(i,j) + prior(i,j); 
    end
end

% poolobj = parpool(n_pools);
% 
% parfor j = 1:n_pools
%     for i = 1:length(queue)
%         ss(i,j) = Arm16_SS(chain(queue(i),:,j),data);
%         prior(i,j) = prior_act(chain(queue(i),:,j),mu,sig);    
%         sschain(i,j) = ss(i,j) + prior(i,j); 
%     end
% end
% delete(poolobj)

toc

output = cat(3,ss,prior,sschain); 

filename = 'sschain_result.mat'; 
save(filename,'output');

color_scheme = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};

%%
figure(1)

set(gcf,'units','centimeters','Position',[7.5935 4.2863 18 18])

subplot(3,1,1)
for j = 1:n_pools
    plot(ss(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('ss')
ylim([0 2500])
ylabel('sum of square error')

subplot(3,1,2)
for j = 1:n_pools
    plot(prior(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('prior')
ylim([0 40])
ylabel('sum musc exc cubed')
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','orientation','horizontal')
legend('boxoff')

subplot(3,1,3)
for j = 1:n_pools
    plot(sschain(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('ss + prior')
ylim([0 2500])
xlabel('iteration')

%
figure(2)

set(gcf,'units','centimeters','Position',[30 4.2863 14 18])

subplot(3,1,1)
for j = 1:n_pools
    plot(ss(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('Likelihood')
text(-0.1,1.1,'A','fontsize',11,'fontweight','bold','units','normalized')
ylim([0 200])
yticks([0 100 200])
ylabel('Sum of squared error')
box off

subplot(3,1,2)
for j = 1:n_pools
    plot(prior(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('Prior')
ylim([0 40])
text(-0.1,1.1,'B','fontsize',11,'fontweight','bold','units','normalized')
ylabel('Sum musc. exc. cubed')
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','orientation','horizontal','numColumns',4)
legend('boxoff')
ylim([0 200])
yticks([0 100 200])
box off

subplot(3,1,3)
for j = 1:n_pools
    plot(sschain(:,j),'color', color_scheme{j},'LineWidth',2)
    hold on 
end
title('Posterior')
text(-0.1,1.1,'C','fontsize',11,'fontweight','bold','units','normalized')
ylim([0 200])
yticks([0 100 200])
xlabel('iteration')
box off

%%

function ss = Arm16_SS(theta,data)
 
%     time   = data.xdata;
%     ydata  = data.ydata(:,2:end);
    ydata_pos = data.ydata(:,2); % this is for only caring about position
    ydata_vel = data.ydata(:,3);
    % we know initial state is pos = 1, vel = 0, but could add this as
    % variability later... 
       
    y = Arm16_fun(theta(1:60));
%     pre_temp = theta(61); 
    
%     ymodel = y(:,3:4);
    ymodel_pos = y(:,1); % position
    ymodel_vel = y(:,2); % velocity
%     activations = y(:,5:2:15);
%     act_sq = activations.^2;
%     act_sq_sum = sum(sum(act_sq));
% %     
    likelihood(1) = sum(((ymodel_pos - ydata_pos).^2)/0.4^2);
    likelihood(2) = sum(((ymodel_vel - ydata_vel).^2)/1.6^2);
        % add in "periodicity" terms to help get it better. 
%     likelihood(2) = (1/3) * ((ymodel_pos(1) - ymodel_pos(end))^2)/0.03^2; 
%     likelihood(3) = (1/3) * ((ymodel_vel(1) - ymodel_vel(end))^2)/0.3^2; 
%     likelihood = (1/2) * (sum(((ymodel - ydata).^2)));
    
    % calcualte the prior here instead
%     tfinal = 0.4;
%     h = 0.001;
%     time = (0:h:tfinal)';
%     Muscle_1 = theta(1:8);
%     Muscle_2 = theta(9:16);
%     Muscle_3 = theta(17:24);
%     Muscle_4 = theta(25:32);
%     Muscle_5 = theta(33:40);
%     Muscle_6 = theta(41:48);
%     
%     controls(:,1) = CRBF_excit(time,Muscle_1);
%     controls(:,2) = CRBF_excit(time,Muscle_2);
%     controls(:,3) = CRBF_excit(time,Muscle_3);
%     controls(:,4) = CRBF_excit(time,Muscle_4);
%     controls(:,5) = CRBF_excit(time,Muscle_5);
%     controls(:,6) = CRBF_excit(time,Muscle_6);
%     
%     SumIntegExc = 0;
%     
%     for k = 1:6
%         SumIntegExc = SumIntegExc + (h * trapz(controls(:,k).^3));
%     end
% 
% %     prior = (SumIntegExc/.0005)^2;
%     prior = ((SumIntegExc - 0.0009)/.0001)^2;
% 
%     prior = (1/2) * (SumIntegExc/.1)^2;


%       offset = 4; % try between 2 and 4
%       amplitude = 0.4; % try between 0 and 1
%     
%     ss = likelihood + prior; 
%     denominator = temperature(pre_temp,offset,amplitude);
    ss = sum(likelihood);  

%      ss = sum(likelihood); 
%     ss = sum(((ymodel - ydata).^2)/0.25^2);


    % convert from radians to degrees, to see if that's better? 
%     ymodel_deg = rad2deg(ymodel);
%     ydata_deg  = rad2deg(ydata);
    
%     ss = sum(((ymodel_deg - ydata_deg).^2)/500);
%     ss = sum(((ymodel - ydata).^2)/0.25^2);

end

% function prior = zero_prior(theta,mu,sig)
%     prior =0; 
% end 

function prior = prior_act(theta,mu,sig)
    tfinal = 0.5;
    h = 0.001;
    time = (0:h:tfinal)';
    Muscle_1 = theta(1:10);
    Muscle_2 = theta(11:20);
    Muscle_3 = theta(21:30);
    Muscle_4 = theta(31:40);
    Muscle_5 = theta(41:50);
    Muscle_6 = theta(51:60);
    
%     pre_temp = theta(61); 
    
    controls(:,1) = CRBF_excit(time,Muscle_1);
    controls(:,2) = CRBF_excit(time,Muscle_2);
    controls(:,3) = CRBF_excit(time,Muscle_3);
    controls(:,4) = CRBF_excit(time,Muscle_4);
    controls(:,5) = CRBF_excit(time,Muscle_5);
    controls(:,6) = CRBF_excit(time,Muscle_6);
    
    SumIntegExc = 0;
    
    for k = 1:6
        SumIntegExc = SumIntegExc + (h * trapz(controls(:,k).^3));
    end

%     prior = (SumIntegExc/.04)^2; 

    prior = ((SumIntegExc - 0.0)/0.0671)^2;
     
    
%     prior2 = (pre_temp)^2;
%     
%     
%     offset = 4; % try between 2 and 4
%     amplitude = 0.4; % try between 0 and 1
%     denominator = temperature(pre_temp,offset,amplitude);
%     
%     prior_sum = (prior1/denominator) + prior2; 
    
    % post calc when temp is between less than 1.001. 
end 

% function data = temperature(pre_temp,offset,amplitude) 
%     
%     data = 1 + inverse_logit(pre_temp-offset) * amplitude; 
%     
% end
% 
% function output = inverse_logit(input)
%     if input>0 
%         output = 1./(1+exp(-input));
%     else
%         output = exp(input)./(1+exp(input));
%     end
% end

function y=Arm16_fun(theta)
% time_100 = [0:0.004:0.4];
%     raw_out = Arm16_forward_controls_CRBF(theta); 
%     y =  interp1(raw_out(:,2),raw_out(:,3),time_100);
    time_100 = [0:0.005:0.5];
    raw_out = Arm16_SimManager_controls_CRBF_6musc(theta); 
    y(:,1) =  interp1(raw_out(:,1),raw_out(:,2),time_100);
    y(:,2) =  interp1(raw_out(:,1),raw_out(:,3),time_100);
    % Need to delete the file because sometimes it doesn't over-write it
    % correctly and sends an error. 
%     delete('C:\Users\rtjoh\Google Drive\Sensitivity_Analyses_Github\Sensitivity-Analysis-MCMC\Arm16_MCMC\Forward_Setup_MCMC\Forward_Setup.xml')
%     trial = trial+1;
    
end