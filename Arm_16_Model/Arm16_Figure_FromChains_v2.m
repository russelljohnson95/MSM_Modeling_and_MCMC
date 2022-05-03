% Plot figures for arm 16 MCMC result
close all
clear
clc


% load('Reference_Data_Arm16_Rigid2.mat')
% 
% position     = ref_data_rigid.positions;
% velocity     = ref_data_rigid.velocities;

% Load in the Static Optimization files (states and forces)
final_time = 0.5; 
time_interval2 = final_time/100;
time_int = (0:time_interval2:final_time)';

% Input_raw = dlmread('arm16_pert4_states.sto','\t',7,0);
input_struct = load('arm16_noise_kinematics.mat');
input_states_ref   = readmatrix('arm16_Tracking_p50_degroote_w75_cubed_v6_states.sto','FileType','text'); 
input_controls_ref = readmatrix('arm16_Tracking_p50_degroote_w75_cubed_v6controls.sto','FileType','text');
input_data = input_struct.new_noise_data; 

position = interp1(input_data(:,1),input_data(:,2),time_int);
velocity = interp1(input_data(:,1),input_data(:,3),time_int);

reference_forces_input = importdata('arm16_forcereport_ForceReporter_forces.sto'); 
reference_forces = reference_forces_input.data; 

% StatOpStates = readmatrix('arm16_pert4_states.sto','FileType','text');
% StatOpForces = readmatrix('arm16_pert4_ForceReporter_forces.sto','FileType','text');

% return

rng(99)


% % Load in the Computed Muscle Controls files (states and forces)
% CMCStates   = readmatrix('arm16_CMC_v3_states.sto','FileType','text');
% CMCForces   = readmatrix('arm16_CMC_v3_Actuation_force.sto','FileType','text');
% 
% % create new time column that takes off the offset that we had to put in
% % for the CMC algorithm
% Interval = 2/300;
% Time = (0:Interval:0.500)';
% 
% TimeCMC = 0.03:Interval:0.53;
% 
% % Need to downsample and interpolate the CMC data to fit the static optimization
% CMCStates_Intp = interp1(CMCStates(1:end-1,1),CMCStates(1:end-1,2:end),TimeCMC);
% CMCForces_Intp = interp1(CMCStates(1:end-1,1),CMCForces(1:end-1,2:end),TimeCMC);


% time_SO  = StatOpStates(:,1); 
% position = interp1(StatOpStates(:,1),StatOpStates(:,2),time_int);
% velocity = interp1(StatOpStates(:,1),StatOpStates(:,3),time_int);

% load chain_results_20201117T012029 %% RUN ONE!!
% load chain_results_20201118T062322 % Narrow Prior Function
% load chain_results_20201119T001613 % intermediate Prior Function
% load chain_results_20201119T225307 
% load chain_results_20201116T045918
% load chain_results_20201212T072650 %% RUN TWO!!
% load chain_results_20210204T235920
% load  chain_results_20210304T114654
load chain_results_20220115T023505
load results_20220115T023505





n_iter = size(chain,1);
n_pools = size(chain,3); 
burn_in = size(chain,1) *0.50;
options.nsimu = size(chain,1);

% chain_corr = chain(:,:,1); 

% for i = 1:60
%     AutoCorr_msd(:,i) = autocorr(chain_corr(:,i),n_iter-1); 
%     AutoCorr_temp = AutoCorr_msd( 1:find( AutoCorr_msd(:,i) < 0, 1 ) );
%     Neff(i) = (length(AutoCorr_temp)*n_iter)/(-1+2*sum(AutoCorr_temp));
% end

% for i = 1:60 % number of parameters = 60
% %     AutoCorr_msd(:,i) = autocorr(chain_corr(:,i)); 
%     AutoCorr_msd(:,i) = autocorr(chain_corr(:,i),n_iter-1); %auto corr of 1 chain, with lags from 1:n-1
%     AutoCorr_temp = AutoCorr_msd( 1:find( AutoCorr_msd(:,i) < 0, 1 ) ); % find the correlation until value goes < 1
% %     Neff(i) = (length(AutoCorr_temp)*n_iter)/(-1+2*sum(AutoCorr_temp));
%     Neff(i) = (length(AutoCorr_temp)*n_iter)/(-1+2*sum(AutoCorr_temp)); % (M * N)/t
% end
% 
% for i = 1:60
%     [mESS(i),Sigma(i),b(i)] = multiESS(chain_corr(:,i));
% end

% reduce the size of the chain to pick out every 100 iterations after
% burn-in

for i = 1:n_pools 
   chain_reduce(:,:,i) = chain(burn_in+1:100:end,:,i);
end



choose = 4; 
names = {'TLo1','TLo2','TLo3','TLo4','TLo5','TLo6','TLo7','TLo8','TLo9','TLo10',...
         'TLa1','TLa2','TLa3','TLa4','TLa5','TLa6','TLa7','TLa8','TLa9','TLa10',...
         'TMe1','TMe2','TMe3','TMe4','TMe5','TMe6','TMe7','TMe8','TMe9','TMe10',...
         'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10',...
         'BS1','BS2','BS3','BS4','BS5','BS6','BS7','BS8','BS9','BS10',...
         'Br1','Br2','Br3','Br4','Br5','Br6','Br7','Br8','Br9','Br10'};

figure(101)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain_reduce(:,1:10,choose),'panellims',names(1:10),[],0)
sgtitle('Tri Long')

figure(102)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain_reduce(:,51:60,choose),'panellims',names(51:60),[],0)
sgtitle('Brachialis')

figure(111)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain_reduce(:,1:10:60,choose),'panellims',names(1:10:60),[],0)
sgtitle('Node 1')


% return

figure(1)

set(gcf,'units','centimeters','Position',[7.5935 4.2863 9 9])

subplot(2,1,1)
% plot(time,position);
% hold on 
plot(time_int,rad2deg(position),'r--','LineWidth',2)
ylabel('position (deg)')
ylim([0 100])
yticks([0 50 100])
xlim([0 0.5])
xticks([0 0.25 0.50])
set(gca,'FontSize',10)
% legend('Reference')
% legend('boxoff')
box off

subplot(2,1,2)
plot(time_int,rad2deg(velocity),'r--','LineWidth',2)
ylabel('velocity (deg/s)')
xlabel('time (s)')
ylim([-800 800])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-800 0 800])
set(gca,'FontSize',10)
box off


addit = 0:0.002:0.2;
% 
position_example = position + addit';

figure(91)

% plot(time,position);
% hold on 
plot(time_int,rad2deg(position),'r--','LineWidth',2)
hold on 
plot(time_int,rad2deg(position_example),'k','LineWidth',2)
% ylabel('position (deg)')
% xlabel('time (s)')
ylim([0 100])
yticks([0 50 100])
set(gca,'FontSize',15)
xlim([0 0.5])
xticks([0 0.25 0.5])
legend('Reference','MCMC')
legend('boxoff')
box off

% return

% choose = 4; % chose chain/result 1

% chain_plot = chain_reduce(:,:,choose);

% load('theta_input3.mat')
% position     = Ref_Data.Position;
% velocity     = Ref_Data.Velocity;
% theta_ref    = Ref_Data.Theta;
% theta_ref = theta_input3;
% n_pools = 7;
% n_draws = 25;

draw = 1:100:size(chain_reduce,1);
n_draws = length(draw); 

for j = 1:n_pools
    for k = 1:n_draws
%        draw(k) = randi([burn_in options.nsimu]);
       Draw_Results(k,:,j) = chain_reduce(draw(k),:,j);
       yy = Arm16_SimManager_controls_CRBF_6musc(Draw_Results(k,:,j));
       positions_draw(:,k,j) = interp1(yy(:,1),yy(:,2),time_int);
       velocities_draw(:,k,j) = interp1(yy(:,1),yy(:,3),time_int);
    end
end

% prepare data for plotting: mean/sd
positions_data = horzcat(positions_draw(:,:,1),positions_draw(:,:,2),positions_draw(:,:,3),positions_draw(:,:,4),positions_draw(:,:,5),positions_draw(:,:,6),positions_draw(:,:,7));
positions_mean = mean(positions_data,2);
positions_std  = std(positions_data,0,[2]);

velocities_data= horzcat(velocities_draw(:,:,1),velocities_draw(:,:,2),velocities_draw(:,:,3),velocities_draw(:,:,4),velocities_draw(:,:,5),velocities_draw(:,:,6),velocities_draw(:,:,7));
velocities_mean = mean(velocities_data,2);
velocities_std  = std(velocities_data,0,[2]);


% calculate RMS error between reference and MCMC results
for i = 1:(n_draws*n_pools)
    RMSE_Position(i) = sqrt(mean((positions_data(:,i)-position).^2));
    RMSE_Velocity(i) = sqrt(mean((velocities_data(:,i)-velocity).^2));
end

RMSE_Position_Mean = mean(RMSE_Position);
RMSE_Velocity_Mean = mean(RMSE_Velocity); 

% return

% 
for j = 1:n_pools
    for i = 1:10 % this 10 is for the number of nodes per muscle. not the number of random draws.. 
        Muscle_draw1(:,i,j) = Draw_Results(:,i,j);
        Muscle_draw2(:,i,j) = Draw_Results(:,i+10,j);
        Muscle_draw3(:,i,j) = Draw_Results(:,i+20,j);
        Muscle_draw4(:,i,j) = Draw_Results(:,i+30,j);
        Muscle_draw5(:,i,j) = Draw_Results(:,i+40,j);
        Muscle_draw6(:,i,j) = Draw_Results(:,i+50,j);
    end
end
tfinal = 0.5;
time = (0:0.001:tfinal)';
for j = 1:n_pools
    for i = 1:n_draws
        Muscle_draw1_traj(:,i,j) = CRBF_excit(time,Muscle_draw1(i,:,j));
        Muscle_draw2_traj(:,i,j) = CRBF_excit(time,Muscle_draw2(i,:,j));
        Muscle_draw3_traj(:,i,j) = CRBF_excit(time,Muscle_draw3(i,:,j));    
        Muscle_draw4_traj(:,i,j) = CRBF_excit(time,Muscle_draw4(i,:,j));
        Muscle_draw5_traj(:,i,j) = CRBF_excit(time,Muscle_draw5(i,:,j));
        Muscle_draw6_traj(:,i,j) = CRBF_excit(time,Muscle_draw6(i,:,j));
    end
end
for j = 1:n_pools
    for i = 1:n_draws
        controls(:,1,i,j) = CRBF_excit(time,Muscle_draw1(i,:,j));
        controls(:,2,i,j) = CRBF_excit(time,Muscle_draw2(i,:,j));
        controls(:,3,i,j) = CRBF_excit(time,Muscle_draw3(i,:,j));    
        controls(:,4,i,j) = CRBF_excit(time,Muscle_draw4(i,:,j));
        controls(:,5,i,j) = CRBF_excit(time,Muscle_draw5(i,:,j));
        controls(:,6,i,j) = CRBF_excit(time,Muscle_draw6(i,:,j));
    end
end

SumIntegExcOut = zeros(n_draws,n_pools);
h = 0.001;

for j = 1:n_pools
    for i = 1:n_draws
        for k = 1:6
            SumIntegExcOut(i,j) = SumIntegExcOut(i,j) + (h * trapz(controls(:,k,i,j).^3));
        end
    end
end

% Create the rank plots to check for convergence 
% rank_plot_Arm16_10CRBVs_fxn(chain(burn_in:end,:,:))
% muscles = {'Tri Long', 'Tri Lat','Tri Med','Biceps LH', 'Biceps SH','Brachior'};
% muscle 4 (Biceps LH)
nNodes = 10;
Chain_rank_plot_muscle_10node(chain_reduce(:,((1*nNodes)-9):(1*nNodes),:),'Triceps Long',burn_in)
Chain_rank_plot_muscle_10node(chain_reduce(:,((2*nNodes)-9):(2*nNodes),:),'Triceps Lat',burn_in)
Chain_rank_plot_muscle_10node(chain_reduce(:,((3*nNodes)-9):(3*nNodes),:),'Triceps Med',burn_in)
Chain_rank_plot_muscle_10node(chain_reduce(:,((4*nNodes)-9):(4*nNodes),:),'Biceps LH',burn_in)
Chain_rank_plot_muscle_10node(chain_reduce(:,((5*nNodes)-9):(5*nNodes),:),'Biceps SH',burn_in)
Chain_rank_plot_muscle_10node(chain_reduce(:,((6*nNodes)-9):(6*nNodes),:),'Brachialis',burn_in)

%%

figure()

color_scheme = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};

set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 8])

histogram(SumIntegExcOut(:,1),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{1},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,2),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{2},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,3),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{3},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,4),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{4},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,5),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{5},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,6),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{6},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,7),'NumBins',40,'BinLimits',[0 .30],'FaceColor',color_scheme{7},'FaceAlpha',0.5);
xlabel('Sum Muscle Excitation Cubed')
set(gca,'FontSize',10)
ylim([0 8])
yticks([0 4 8])
box off
xlim([0 0.4])
xticks([0 0.1 0.2 0.3 0.4])
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','location','northeast')
legend('boxoff')

SumIntegExcOut_Vector = vertcat(SumIntegExcOut(:,1),SumIntegExcOut(:,2),SumIntegExcOut(:,3),SumIntegExcOut(:,4),SumIntegExcOut(:,5),SumIntegExcOut(:,6),SumIntegExcOut(:,7));


figure()

set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 8])

pd1 = histfit(SumIntegExcOut(:,1),[],'normal');
pd1(1).FaceColor = 'none';
pd1(1).EdgeColor = 'none';
pd1(2).Color     = color_scheme{1};
hold on 
pd2 = histfit(SumIntegExcOut(:,2),[],'normal');
pd2(1).FaceColor = 'none';
pd2(1).EdgeColor = 'none';
pd2(2).Color     = color_scheme{2};
hold on 
pd3 = histfit(SumIntegExcOut(:,3),[],'normal');
pd3(1).FaceColor = 'none';
pd3(1).EdgeColor = 'none';
pd3(2).Color     = color_scheme{3};
hold on 
pd4 = histfit(SumIntegExcOut(:,4),[],'normal');
pd4(1).FaceColor = 'none';
pd4(1).EdgeColor = 'none';
pd4(2).Color     = color_scheme{4};
hold on 
pd5 = histfit(SumIntegExcOut(:,5),[],'normal');
pd5(1).FaceColor = 'none';
pd5(1).EdgeColor = 'none';
pd5(2).Color     = color_scheme{5};
hold on 
pd6 = histfit(SumIntegExcOut(:,6),[],'normal');
pd6(1).FaceColor = 'none';
pd6(1).EdgeColor = 'none';
pd6(2).Color     = color_scheme{6};
hold on 
pd7 = histfit(SumIntegExcOut(:,7),[],'normal');
pd7(1).FaceColor = 'none';
pd7(1).EdgeColor = 'none';
pd7(2).Color     = color_scheme{7};
xlim([0 0.4])
xlabel('Sum Muscle Excitation Cubed')
set(gca,'FontSize',10)
% ylim([0 8])
% yticks([0 4 8])
box off
xticks([0 0.1 0.2 0.3 0.4])
legend([pd1(2),pd2(2),pd3(2),pd4(2),pd5(2),pd6(2),pd7(2)],'Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','location','northeast')
legend('boxoff')

% SumIntegExcOut_Vector = vertcat(SumIntegExcOut(:,1),SumIntegExcOut(:,2),SumIntegExcOut(:,3),SumIntegExcOut(:,4));

%% correlations

% Plot the correlations within a muscle 
names = {'TLo1','TLo2','TLo3','TLo4','TLo5','TLo6','TLo7','TLo8','TLo9','TLo10',...
         'TLa1','TLa2','TLa3','TLa4','TLa5','TLa6','TLa7','TLa8','TLa9','TLa10',...
         'TMe1','TMe2','TMe3','TMe4','TMe5','TMe6','TMe7','TMe8','TMe9','TMe10',...
         'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10',...
         'BS1','BS2','BS3','BS4','BS5','BS6','BS7','BS8','BS9','BS10',...
         'Br1','Br2','Br3','Br4','Br5','Br6','Br7','Br8','Br9','Br10'};

figure(101)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain_reduce(:,1:10,choose),'panellims',names(1:10),[],0)
sgtitle('Tri Long')
% 
% figure(102)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,11:20,choose),'panellims',names(11:20),[],0)
% sgtitle('Tri Lat')
% 
% figure(103)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,21:30,choose),'panellims',names(21:30),[],0)
% sgtitle('Tri Med')
% 
% figure(104)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,31:40,choose),'panellims',names(31:40),[],0)
% sgtitle('Biceps LH')
% 
% figure(105)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,41:50,choose),'panellims',names(41:50),[],0)
% sgtitle('Biceps SH')
% 
% figure(106)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,51:60,choose),'panellims',names(51:60),[],0)
% sgtitle('Brachioradialis')

%  Plot the correlations within a node
figure(111)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain_reduce(:,1:10:60,choose),'panellims',names(1:10:60),[],0)
sgtitle('Node 1')

% figure(112)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,2:10:60,choose),'panellims',names(2:10:60),[],0)
% sgtitle('Node 2')
% 
% figure(113)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,3:10:60,choose),'panellims',names(3:10:60),[],0)
% sgtitle('Node 3')
% 
% figure(114)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,4:10:60,choose),'panellims',names(4:10:60),[],0)
% sgtitle('Node 4')
% 
% figure(115)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,5:10:60,choose),'panellims',names(5:10:60),[],0)
% sgtitle('Node 5')
% 
% figure(116)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,6:10:60,choose),'panellims',names(6:10:60),[],0)
% sgtitle('Node 6')
% 
% figure(117)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,7:10:60,choose),'panellims',names(7:10:60),[],0)
% sgtitle('Node 7')
% 
% figure(118)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,8:10:60,choose),'panellims',names(8:10:60),[],0)
% sgtitle('Node 8')
% 
% figure(119)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain_reduce(burn_in:end,9:10:60,choose),'panellims',names(9:10:60),[],0)
% sgtitle('Node 9')
% 
% figure(120)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
% pairs_v3(chain_reduce(burn_in:end,10:10:60,choose),'panellims',names(10:10:60),[],0)
% sgtitle('Node 10')



%%

% return

figure()
pd8 = histfit(SumIntegExcOut_Vector,[],'beta');
pd8(1).FaceColor = 'none';
pd8(1).EdgeColor = 'none';
pd8(2).Color     = [.2 .2 .2];
xlim([0 0.3])
h2 = xline(0.07,'--r','LineWidth',1.5);
xlabel('Sum Muscle Exc. Cubed')
ylabel('Density')
set(gca,'FontSize',15)
legend([pd8(2),h2],'density','target')
set(gca,'fontsize',16)

% pull out muscle excitations from file
% muscle_1_controls = CRBF_excit(time,theta_ref(1:10));
% muscle_2_controls = CRBF_excit(time,theta_ref(11:20));
% muscle_3_controls = CRBF_excit(time,theta_ref(21:30));
% muscle_4_controls = CRBF_excit(time,theta_ref(31:40));
% muscle_5_controls = CRBF_excit(time,theta_ref(41:50));
% muscle_6_controls = CRBF_excit(time,theta_ref(51:60));


chain_input = vertcat(Draw_Results(:,:,1),Draw_Results(:,:,2),Draw_Results(:,:,3),Draw_Results(:,:,4),Draw_Results(:,:,5),Draw_Results(:,:,6),Draw_Results(:,:,7));

for i = 1:size(chain_input,1)
    [kinematics_out(:,:,i),force_out(:,:,i)] = Arm16_SimManager_controls_CRBF_6musc_wForce(chain_input(i,:));
end

Muscle_1_traj_mean = mean(force_out(:,1,:),3);
Muscle_2_traj_mean = mean(force_out(:,2,:),3);
Muscle_3_traj_mean = mean(force_out(:,3,:),3);
Muscle_4_traj_mean = mean(force_out(:,4,:),3);
Muscle_5_traj_mean = mean(force_out(:,5,:),3);
Muscle_6_traj_mean = mean(force_out(:,6,:),3);

Muscle_1_traj_std(:,1) = std(force_out(:,1,:),0,[3]);
Muscle_2_traj_std(:,1) = std(force_out(:,2,:),0,[3]);
Muscle_3_traj_std(:,1) = std(force_out(:,3,:),0,[3]);
Muscle_4_traj_std(:,1) = std(force_out(:,4,:),0,[3]);
Muscle_5_traj_std(:,1) = std(force_out(:,5,:),0,[3]);
Muscle_6_traj_std(:,1) = std(force_out(:,6,:),0,[3]);

% figure(101)
% 
% subplot(3,2,1)
% [Pa,Li] = JackKnife(time_int,Muscle_1_traj_mean,Muscle_1_traj_std,'k',[0.75, 0.75, 0.75]);
% ylim([0 600])

% return

% Muscle_draw1_traj_mean = mean(mean(Muscle_draw1_traj,2),3);
% Muscle_draw2_traj_mean = mean(mean(Muscle_draw2_traj,2),3);
% Muscle_draw3_traj_mean = mean(mean(Muscle_draw3_traj,2),3);
% Muscle_draw4_traj_mean = mean(mean(Muscle_draw4_traj,2),3);
% Muscle_draw5_traj_mean = mean(mean(Muscle_draw5_traj,2),3);
% Muscle_draw6_traj_mean = mean(mean(Muscle_draw6_traj,2),3);
% 
% Muscle_draw1_traj_std = std(Muscle_draw1_traj,0,[2,3]);
% Muscle_draw2_traj_std = std(Muscle_draw2_traj,0,[2,3]);
% Muscle_draw3_traj_std = std(Muscle_draw3_traj,0,[2,3]);
% Muscle_draw4_traj_std = std(Muscle_draw4_traj,0,[2,3]);
% Muscle_draw5_traj_std = std(Muscle_draw5_traj,0,[2,3]);
% Muscle_draw6_traj_std = std(Muscle_draw6_traj,0,[2,3]);
% 
% Muscle_draw1_traj_std_min = Muscle_draw1_traj_mean - Muscle_draw1_traj_std;
% Muscle_draw2_traj_std_min = Muscle_draw2_traj_mean - Muscle_draw2_traj_std;
% Muscle_draw3_traj_std_min = Muscle_draw3_traj_mean - Muscle_draw3_traj_std;
% Muscle_draw4_traj_std_min = Muscle_draw4_traj_mean - Muscle_draw4_traj_std;
% Muscle_draw5_traj_std_min = Muscle_draw5_traj_mean - Muscle_draw5_traj_std;
% Muscle_draw6_traj_std_min = Muscle_draw6_traj_mean - Muscle_draw6_traj_std;
% 
% Muscle_draw1_traj_std_min ( Muscle_draw1_traj_std_min < 0 ) = 0;
% Muscle_draw2_traj_std_min ( Muscle_draw2_traj_std_min < 0 ) = 0;
% Muscle_draw3_traj_std_min ( Muscle_draw3_traj_std_min < 0 ) = 0;
% Muscle_draw4_traj_std_min ( Muscle_draw4_traj_std_min < 0 ) = 0;
% Muscle_draw5_traj_std_min ( Muscle_draw5_traj_std_min < 0 ) = 0;
% Muscle_draw6_traj_std_min ( Muscle_draw6_traj_std_min < 0 ) = 0;
% 
% Muscle_draw1_traj_std_max = Muscle_draw1_traj_mean + Muscle_draw1_traj_std;
% Muscle_draw2_traj_std_max = Muscle_draw2_traj_mean + Muscle_draw2_traj_std;
% Muscle_draw3_traj_std_max = Muscle_draw3_traj_mean + Muscle_draw3_traj_std;
% Muscle_draw4_traj_std_max = Muscle_draw4_traj_mean + Muscle_draw4_traj_std;
% Muscle_draw5_traj_std_max = Muscle_draw5_traj_mean + Muscle_draw5_traj_std;
% Muscle_draw6_traj_std_max = Muscle_draw6_traj_mean + Muscle_draw6_traj_std;



% so take the draw results and concatenate them vertically
% chain_fromDraw = vertcat(Draw_Results(:,:,1),Draw_Results(:,:,2),Draw_Results(:,:,3),Draw_Results(:,:,4));
for j = 1:n_pools
    for i = 1:n_draws
        controls_input(:,(j*n_draws-(n_draws)+i),1) = CRBF_excit(time_int,Muscle_draw1(i,:,j));
        controls_input(:,(j*n_draws-(n_draws)+i),2) = CRBF_excit(time_int,Muscle_draw2(i,:,j));
        controls_input(:,(j*n_draws-(n_draws)+i),3) = CRBF_excit(time_int,Muscle_draw3(i,:,j));    
        controls_input(:,(j*n_draws-(n_draws)+i),4) = CRBF_excit(time_int,Muscle_draw4(i,:,j));
        controls_input(:,(j*n_draws-(n_draws)+i),5) = CRBF_excit(time_int,Muscle_draw5(i,:,j));
        controls_input(:,(j*n_draws-(n_draws)+i),6) = CRBF_excit(time_int,Muscle_draw6(i,:,j));
    end
end

% arm16_contourplot_6mus(chain_fromDraw,time,time_int,[],[],[],[],[])
% arm16_contoutplot_6mus_v2(controls_input,time_int);

% figure()
% 
% subplot(2,3,1)
% [Pa,Li] = JackKnife(time_int,Muscle_1_traj_mean,Muscle_1_traj_std,'k',[0.75, 0.75, 0.75]);
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,1,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,1,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,1,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,1,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% legend([h1 h2 h3 h4],'Chain 1','Chain 2','Chain 3','Chain 4','Orientation','horizontal')
% legend('boxoff')
% box off
% set(gca,'fontsize',16)
% title('Triceps Long.')
% ylabel('Force (N)')
% ylim([0 800])
% xlim([0 0.5])
% 
% % -----
% subplot(2,3,2)
% [Pa,Li] = JackKnife(time_int,Muscle_2_traj_mean,Muscle_2_traj_std,'k',[0.75, 0.75, 0.75]);
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,2,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,2,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,2,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,2,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% box off
% set(gca,'fontsize',16)
% title('Triceps Lat.')
% % ylabel('Force (N)')
% ylim([0 625])
% xlim([0 0.5])
% 
% subplot(2,3,3)
% [Pa,Li] = JackKnife(time_int,Muscle_3_traj_mean,Muscle_3_traj_std,'k',[0.75, 0.75, 0.75]);
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,3,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,3,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,3,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,3,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% ylim([0 625])
% xlim([0 0.5])
% box off
% set(gca,'fontsize',16)
% title('Triceps Med.')
% 
% 
% subplot(2,3,4)
% [Pa,Li] = JackKnife(time_int,Muscle_4_traj_mean,Muscle_4_traj_std,'k',[0.75, 0.75, 0.75]);
% 
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,4,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,4,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,4,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,4,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% ylim([0 625])
% xlim([0 0.5])
% box off
% set(gca,'fontsize',16)
% title('Biceps LH')
% ylabel('Force (N)')
% xlabel('time (s)')
% 
% 
% 
% subplot(2,3,5)
% [Pa,Li] = JackKnife(time_int,Muscle_5_traj_mean,Muscle_5_traj_std,'k',[0.75, 0.75, 0.75]);
% 
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,5,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,5,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,5,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,5,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% ylim([0 500])
% xlim([0 0.5])
% box off
% set(gca,'fontsize',16)
% title('Biceps SH')
% % ylabel('Force (N)')
% xlabel('time (s)')
% 
% subplot(2,3,6)
% [Pa,Li] = JackKnife(time_int,Muscle_6_traj_mean,Muscle_6_traj_std,'k',[0.75, 0.75, 0.75]);
% 
% hold on 
% for i = 1:n_draws
%     h4 = plot(time_int,force_out(:,6,i+3*n_draws),'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.0);
%     h4.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h1 = plot(time_int,force_out(:,6,i),'color',[0, 0.5, 0],'LineWidth',1.0);
%     h1.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h2 = plot(time_int,force_out(:,6,i+n_draws,1),'color',[0.75, 0.75, 0],'LineWidth',1.0);
%     h2.Color(4) = 0.15;
%     hold on 
% end
% for i = 1:n_draws
%     h3 = plot(time_int,force_out(:,6,i+2*n_draws,1),'color',[0, 0.75, 0.75],'LineWidth',1.0);
%     h3.Color(4) = 0.15;
%     hold on 
% end
% ylim([0 1100])
% xlim([0 0.5])
% box off
% set(gca,'fontsize',16)
% title('Brachialis')
% % ylabel('Force (N)')
% xlabel('time (s)')


%%


figure(21)

set(gcf,'units','centimeters','Position',[7.5935 4.2863 21 15])

% this section is for the kinematics (positions and velocities)
subplot(3,3,1)
[Pa,Li] = JackKnife(time_int,rad2deg(positions_mean),rad2deg(positions_std),'k',[0.75, 0.75, 0.75]);
hold on 
h2 = plot(time_int,rad2deg(position),'r--','LineWidth',1.5);
set(gca,'fontsize',10)
% h2 = plot(time_int,position,'r--','LineWidth',1.5);
ylabel('position (deg)')
% xlabel('time (s)')
text(-0.3,1.1,'A','fontsize',10,'fontweight','bold','units','normalized')
ylim([0 100])
yticks([0 50 100])
xlim([0 0.5])
xticks([0 0.25 0.50])
% legend([Li,h2],"MCMC","Reference","Location","Northwest",'orientation','horizontal');
% legend('boxoff')
box off



subplot(3,3,2)

[Pa,Li] = JackKnife(time_int,rad2deg(velocities_mean),rad2deg(velocities_std),'k',[0.75, 0.75, 0.75]);
set(gca,'fontsize',10)
h2 = plot(time_int,rad2deg(velocity),'r--','LineWidth',1.5);
% h2 = plot(time_int,velocity,'r--','LineWidth',1.5);
ylabel('velocity (deg/s)')
% xlabel('time (s)')
text(-0.3,1.1,'B','fontsize',10,'fontweight','bold','units','normalized')
ylim([-800 800])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-800 0 800])
% legend([h1,h2],"Random Draw","Reference","Location","Northwest");
box off

% This section is the muscle excitations cubed
% first set up the prior distribution
mu  = 0.00; % center
sig = 0.08; % width
ma = 0.0;
mi = 0.3;
xp = linspace(max(mi,mu-3*sig),min([ma,mu+3*sig]));
yp = norpf(xp,mu,sig^2);
yn = nordf((mi-mu)/sig)+1-nordf((ma-mu)/sig); % area outside bounds

subplot(3,3,3)

h1 = plot(xp,yp,'--','Color','#5C74F5','LineWidth',1.5);
hold on
pd1 = histfit(SumIntegExcOut_Vector,[],'normal');
pd1(1).FaceColor = 'none';
pd1(1).EdgeColor = 'none';
pd1(2).Color     = '#5C74F5';
legend([h1,pd1(2)],'Prior Density','Post. Density');
xlim([0 0.3])
ylim([0 40])
xlabel('Sum Muscle Exc. Cubed')
ylabel('density')
text(-0.3,1.1,'C','fontsize',10,'fontweight','bold','units','normalized')
% legend([pd(2),h2],'density','target')
legend('boxoff')
set(gca,'fontsize',10)
box off



% calculate the width of posterior: 

y_max = max(pd8(2).YData);
B = y_max == pd8(2).YData;
index = find([B] == 1);
posterior_center = pd8(2).XData(index); 

integrate = zeros(100,1);
integrate_sum = zeros(100,1);

for j = index:100
        integrate(j-index+1) = (pd8(2).XData(j)-pd8(2).XData(j-1)) * pd8(2).YData(j);
end

% integrate_sum(1) = integrate(1);
for k = 2:100
    integrate_sum(k) = integrate_sum(k-1)+integrate(k);
end
C = integrate_sum<0.341;
width = index + sum(C);
posterior_std = pd8(2).XData(width)-posterior_center;

% posterior_center
% posterior_std
% 
% return

% the bottom two rows are for the muscle force trajectories. 
subplot(3,3,4)
[Pa,Li] = JackKnife(time_int,Muscle_1_traj_mean,Muscle_1_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(reference_forces(:,1),reference_forces(:,2),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,1)),'b','LineWidth',2)
% xlabel('time (s)')
box off
set(gca,'fontsize',10)
title('Triceps Long.')
text(-0.3,1.1,'D','fontsize',10,'fontweight','bold','units','normalized')
ylabel('force (N)')
ylim([0 800])
yticks([0 400 800])
xlim([0 0.5])
xticks([0 0.25 0.50])

% -----
subplot(3,3,5)
[Pa,Li] = JackKnife(time_int,Muscle_2_traj_mean,Muscle_2_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(reference_forces(:,1),reference_forces(:,3),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,2)),'b','LineWidth',2)
% xlabel('time (s)')
box off
set(gca,'fontsize',10)
title('Triceps Lat.')
text(-0.3,1,'E','fontsize',10,'fontweight','bold','units','normalized')
% ylabel('Force (N)')
ylim([0 625])
yticks([0 300 600])
xlim([0 0.5])
xticks([0 0.25 0.50])


subplot(3,3,6)
[Pa,Li] = JackKnife(time_int,Muscle_3_traj_mean,Muscle_3_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
h1 = plot(reference_forces(:,1),reference_forces(:,4),'--r','LineWidth',2);
% hold on 
% h2 = plot(Time,(CMCForces_Intp(:,3)),'b','LineWidth',2);
% xlabel('time (s)')
ylim([0 625])
yticks([0 300 600])
text(-0.3,1.1,'F','fontsize',10,'fontweight','bold','units','normalized')
box off
set(gca,'fontsize',10)
title('Triceps Med.')
xlim([0 0.5])
xticks([0 0.25 0.50])

subplot(3,3,7)
[Pa,Li] = JackKnife(time_int,Muscle_4_traj_mean,Muscle_4_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(reference_forces(:,1),reference_forces(:,5),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,4)),'b','LineWidth',2)
ylim([0 625])
yticks([0 300 600])
legend([Pa],'');
legend boxoff
box off
set(gca,'fontsize',10)
text(-0.3,1.1,'G','fontsize',10,'fontweight','bold','units','normalized')
title('Biceps LH')
ylabel('force (N)')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])


subplot(3,3,8)
[Pa,Li] = JackKnife(time_int,Muscle_5_traj_mean,Muscle_5_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(reference_forces(:,1),reference_forces(:,6),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,5)),'b','LineWidth',2)
ylim([0 500])
yticks([0 250 500])
legend([Li,h1],'MCMC Mean +/- SD','Reference','orientation','vertical');
% labelhandles(4).FaceColor = [0.75,0.75,0.75];
% labelhandles(3).YData = [0.83 0.83];
% labelhandles(4).XData = [0.0460 0.25]; labelhandles(4).YData = [0.760 0.760];
legend boxoff
box off
set(gca,'fontsize',10)
title('Biceps SH')
text(-0.3,1.1,'H','fontsize',10,'fontweight','bold','units','normalized')
% ylabel('Force (N)')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])

subplot(3,3,9)
[Pa,Li] = JackKnife(time_int,Muscle_6_traj_mean,Muscle_6_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(reference_forces(:,1),reference_forces(:,7),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,6)),'b','LineWidth',2)
text(-0.3,1.1,'I','fontsize',10,'fontweight','bold','units','normalized')
ylim([0 1100])
yticks([0 550 1100])
box off
set(gca,'fontsize',10)
title('Brachialis')
% ylabel('Force (N)')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])



%%
% return

% plot all muscle forces 

figure()


set(gcf,'units','centimeters','Position',[7.5935 4.2863 21 10])

subplot(2,3,1)
% [fillhandle,msg]=jbfill(time,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time,muscle_1_controls,'c','LineWidth',2)
% hold on 
plot(reference_forces(:,1),reference_forces(:,2),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,1,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,1,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,1,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,1,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,1,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,1,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,1,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end

legend([h1 h2 h3 h4 h5 h6 h7],'Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','Orientation','horizontal')
legend('boxoff')
box off
set(gca,'fontsize',10)
title('Triceps Long.')
xticks([0 0.25 0.50])
ylabel('force (N)')
ylim([0 800])
yticks([0 400 800])

subplot(2,3,2)
% [fillhandle,msg]=jbfill(time_int,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time_int,muscle_2_controls,'b','LineWidth',2)
% hold on 
plot(reference_forces(:,1),reference_forces(:,3),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,2,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,2,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,2,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,2,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,2,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,2,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,2,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Triceps Lat.')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 625])
yticks([0 300 600])

subplot(2,3,3)
% [fillhandle,msg]=jbfill(time_int,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time_int,muscle_3_controls,'b','LineWidth',2)
% hold on 
plot(reference_forces(:,1),reference_forces(:,4),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,3,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,3,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,3,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,3,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,3,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,3,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,3,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Triceps Med.')
xticks([0 0.25 0.50])
% xlabel('time_int (s)')
% ylabel('excitation')
ylim([0 625])
yticks([0 300 600])

subplot(2,3,4)
% [fillhandle,msg]=jbfill(time_int,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time_int,muscle_4_controls,'b','LineWidth',2)
% hold on 
plot(reference_forces(:,1),reference_forces(:,5),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,4,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,4,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,4,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,4,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,4,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,4,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,4,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Biceps LH')
ylabel('force (N)')
xticks([0 0.25 0.50])
xlabel('time (s)')
ylim([0 625])
yticks([0 300 600])

subplot(2,3,5)
% [fillhandle,msg]=jbfill(time,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time,muscle_5_controls,'b','LineWidth',2)
% hold on 
plot(reference_forces(:,1),reference_forces(:,6),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,5,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,5,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,5,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,5,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,5,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,5,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,5,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Biceps SH')
xlabel('time (s)')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 500])
yticks([0 250 500])

subplot(2,3,6)
% [fillhandle,msg]=jbfill(time,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time,muscle_6_controls,'b','LineWidth',2)
plot(reference_forces(:,1),reference_forces(:,7),'--r','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time_int,force_out(:,6,i),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time_int,force_out(:,6,i+25),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time_int,force_out(:,6,i+50),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time_int,force_out(:,6,i+75),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time_int,force_out(:,6,i+100),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time_int,force_out(:,6,i+125),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time_int,force_out(:,6,i+150),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Brachialis')
xlabel('time (s)')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 1100])
yticks([0 550 1100])

% plot all muscle excitations

figure()


set(gcf,'units','centimeters','Position',[7.5935 4.2863 21 10])

subplot(2,3,1)
% [fillhandle,msg]=jbfill(time,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time,muscle_1_controls,'c','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw1_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw1_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw1_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw1_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw1_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw1_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw1_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end

legend([h1 h2 h3 h4 h5 h6 h7],'Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','Orientation','horizontal')
legend('boxoff')
box off
set(gca,'fontsize',10)
title('Triceps Long.')
xticks([0 0.25 0.50])
ylabel('excitation')
ylim([0 1])

subplot(2,3,2)
% [fillhandle,msg]=jbfill(time,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time,muscle_2_controls,'b','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw2_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw2_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw2_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw2_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw2_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw2_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw2_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Triceps Lat.')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 1])

subplot(2,3,3)
% [fillhandle,msg]=jbfill(time,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time,muscle_3_controls,'b','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw3_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw3_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw3_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw3_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw3_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw3_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw3_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Triceps Med.')
xticks([0 0.25 0.50])
% xlabel('time (s)')
% ylabel('excitation')
ylim([0 1])


subplot(2,3,4)
% [fillhandle,msg]=jbfill(time,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time,muscle_4_controls,'b','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw4_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw4_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw4_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw4_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw4_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw4_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw4_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Biceps LH')
ylabel('excitation')
xticks([0 0.25 0.50])
xlabel('time (s)')
ylim([0 1])

subplot(2,3,5)
% [fillhandle,msg]=jbfill(time,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time,muscle_5_controls,'b','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw5_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw5_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw5_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw5_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw5_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw5_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw5_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Biceps SH')
xlabel('time (s)')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 1])

subplot(2,3,6)
% [fillhandle,msg]=jbfill(time,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time,muscle_6_controls,'b','LineWidth',2)
hold on 
for i = 1:n_draws
    h4 = plot(time,Muscle_draw6_traj(:,i,4),'color',color_scheme{1},'LineWidth',1.5);
    h4.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h1 = plot(time,Muscle_draw6_traj(:,i,1),'color',color_scheme{2},'LineWidth',1.5);
    h1.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h2 = plot(time,Muscle_draw6_traj(:,i,2),'color',color_scheme{3},'LineWidth',1.5);
    h2.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h3 = plot(time,Muscle_draw6_traj(:,i,3),'color',color_scheme{4},'LineWidth',1.5);
    h3.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h5 = plot(time,Muscle_draw6_traj(:,i,5),'color',color_scheme{5},'LineWidth',1.5);
    h5.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h6 = plot(time,Muscle_draw6_traj(:,i,6),'color',color_scheme{6},'LineWidth',1.5);
    h6.Color(4) = 0.33;
    hold on 
end
for i = 1:n_draws
    h7 = plot(time,Muscle_draw6_traj(:,i,7),'color',color_scheme{7},'LineWidth',1.5);
    h7.Color(4) = 0.33;
    hold on 
end
box off
set(gca,'fontsize',10)
title('Brachialis')
xlabel('time (s)')
xticks([0 0.25 0.50])
% ylabel('excitation')
ylim([0 1])