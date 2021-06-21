% Plot a bunch of results figures for Arm 16 MCMC result

% Author: Russell T. Johnson, rtjohnso@usc.edu
% Last Edited: 6-18-21

% because of the number of random draws and the OpenSim integration,
% running this code takes a few minutes 
% ---------------------------------------------------

close all
clear
clc

% Load in the files from the reference motion 
ReferenceStates = readmatrix('arm16_pert4_states.sto','FileType','text');
ReferenceForces = readmatrix('arm16_pert4_ForceReporter_forces.sto','FileType','text');

rng(99)

final_time = 0.5; 
time_interval2 = final_time/100;
time_int = (0:time_interval2:final_time)';

time_SO  = ReferenceStates(:,1); 
position = interp1(ReferenceStates(:,1),ReferenceStates(:,2),time_int);
velocity = interp1(ReferenceStates(:,1),ReferenceStates(:,3),time_int);

% load results from chain
load  chain_results_20210304T114654

burn_in = size(chain,1) *0.20;
options.nsimu = size(chain,1);

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
ylim([-600 600])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-600 0 600])
set(gca,'FontSize',10)
box off

% return

line = 0:0.002:0.2;

position_example = position + line';

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

choose = 4; % chose chain/result 1

chain_plot = chain(:,:,choose);

load('theta_input3.mat')
% position     = Ref_Data.Position;
% velocity     = Ref_Data.Velocity;
% theta_ref    = Ref_Data.Theta;
theta_ref = theta_input3;
n_pools = 5;
n_draws = 25;

for j = 1:n_pools
    for k = 1:n_draws
       draw(k) = randi([burn_in options.nsimu]);
       Draw_Results(k,:,j) = chain(draw(k),:,j);
       yy = Arm16_SimManager_controls_CRBF_6musc(Draw_Results(k,:,j));
       positions_draw(:,k,j) = interp1(yy(:,1),yy(:,2),time_int);
       velocities_draw(:,k,j) = interp1(yy(:,1),yy(:,3),time_int);
    end
end

% prepare data for plotting: mean/sd
positions_data = horzcat(positions_draw(:,:,1),positions_draw(:,:,2),positions_draw(:,:,3),positions_draw(:,:,4));
positions_mean = mean(positions_data,2);
positions_std  = std(positions_data,0,[2]);

velocities_data= horzcat(velocities_draw(:,:,1),velocities_draw(:,:,2),velocities_draw(:,:,3),velocities_draw(:,:,4));
velocities_mean = mean(velocities_data,2);
velocities_std  = std(velocities_data,0,[2]);

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
Chain_rank_plot_muscle_10node(chain(:,((1*nNodes)-9):(1*nNodes),:),'Triceps Long',burn_in)
Chain_rank_plot_muscle_10node(chain(:,((2*nNodes)-9):(2*nNodes),:),'Triceps Lat',burn_in)
Chain_rank_plot_muscle_10node(chain(:,((3*nNodes)-9):(3*nNodes),:),'Triceps Med',burn_in)
Chain_rank_plot_muscle_10node(chain(:,((4*nNodes)-9):(4*nNodes),:),'Biceps LH',burn_in)
Chain_rank_plot_muscle_10node(chain(:,((5*nNodes)-9):(5*nNodes),:),'Biceps SH',burn_in)
Chain_rank_plot_muscle_10node(chain(:,((6*nNodes)-9):(6*nNodes),:),'Brachialis',burn_in)

%%

figure()

set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 8])

histogram(SumIntegExcOut(:,1),'NumBins',40,'BinLimits',[0 .30],'FaceColor','r','FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,2),'NumBins',40,'BinLimits',[0 .30],'FaceColor','b','FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,3),'NumBins',40,'BinLimits',[0 .30],'FaceColor','g','FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,4),'NumBins',40,'BinLimits',[0 .30],'FaceColor','k','FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,5),'NumBins',40,'BinLimits',[0 .30],'FaceColor','c','FaceAlpha',0.5);
xlabel('Sum Muscle Excitation Cubed')
set(gca,'FontSize',10)
ylim([0 8])
yticks([0 4 8])
box off
xticks([0 0.1 0.2 0.3])
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','location','northeast')
legend('boxoff')

% SumIntegExcOut_Vector = vertcat(SumIntegExcOut(:,1),SumIntegExcOut(:,2),SumIntegExcOut(:,3),SumIntegExcOut(:,4));


figure()

set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 8])

pd1 = histfit(SumIntegExcOut(:,1),[],'normal');
pd1(1).FaceColor = 'none';
pd1(1).EdgeColor = 'none';
pd1(2).Color     = 'r';
hold on 
pd2 = histfit(SumIntegExcOut(:,2),[],'normal');
pd2(1).FaceColor = 'none';
pd2(1).EdgeColor = 'none';
pd2(2).Color     = 'b';
hold on 
pd3 = histfit(SumIntegExcOut(:,3),[],'normal');
pd3(1).FaceColor = 'none';
pd3(1).EdgeColor = 'none';
pd3(2).Color     = 'g';
hold on 
pd4 = histfit(SumIntegExcOut(:,4),[],'normal');
pd4(1).FaceColor = 'none';
pd4(1).EdgeColor = 'none';
pd4(2).Color     = 'k';
hold on 
pd5 = histfit(SumIntegExcOut(:,5),[],'normal');
pd5(1).FaceColor = 'none';
pd5(1).EdgeColor = 'none';
pd5(2).Color     = 'c';
xlim([0 0.3])
xlabel('Sum Muscle Excitation Cubed')
set(gca,'FontSize',10)
% ylim([0 8])
% yticks([0 4 8])
box off
xticks([0 0.1 0.2 0.3])
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','location','northeast')
legend('boxoff')

SumIntegExcOut_Vector = vertcat(SumIntegExcOut(:,1),SumIntegExcOut(:,2),SumIntegExcOut(:,3),SumIntegExcOut(:,4));

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
pairs_v3(chain(burn_in:end,1:10,choose),'panellims',names(1:10),[],0)
sgtitle('Tri Long')
% 
% figure(102)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,11:20,choose),'panellims',names(11:20),[],0)
% sgtitle('Tri Lat')
% 
% figure(103)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,21:30,choose),'panellims',names(21:30),[],0)
% sgtitle('Tri Med')
% 
% figure(104)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,31:40,choose),'panellims',names(31:40),[],0)
% sgtitle('Biceps LH')
% 
% figure(105)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,41:50,choose),'panellims',names(41:50),[],0)
% sgtitle('Biceps SH')
% 
% figure(106)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,51:60,choose),'panellims',names(51:60),[],0)
% sgtitle('Brachioradialis')

%  Plot the correlations within a node
figure(111)
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
pairs_v3(chain(burn_in:end,1:10:60,choose),'panellims',names(1:10:60),[],0)
sgtitle('Node 1')

% figure(112)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,2:10:60,choose),'panellims',names(2:10:60),[],0)
% sgtitle('Node 2')
% 
% figure(113)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,3:10:60,choose),'panellims',names(3:10:60),[],0)
% sgtitle('Node 3')
% 
% figure(114)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,4:10:60,choose),'panellims',names(4:10:60),[],0)
% sgtitle('Node 4')
% 
% figure(115)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,5:10:60,choose),'panellims',names(5:10:60),[],0)
% sgtitle('Node 5')
% 
% figure(116)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,6:10:60,choose),'panellims',names(6:10:60),[],0)
% sgtitle('Node 6')
% 
% figure(117)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,7:10:60,choose),'panellims',names(7:10:60),[],0)
% sgtitle('Node 7')
% 
% figure(118)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,8:10:60,choose),'panellims',names(8:10:60),[],0)
% sgtitle('Node 8')
% 
% figure(119)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 20])
% pairs_v3(chain(burn_in:end,9:10:60,choose),'panellims',names(9:10:60),[],0)
% sgtitle('Node 9')
% 
% figure(120)
% set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 16])
% pairs_v3(chain(burn_in:end,10:10:60,choose),'panellims',names(10:10:60),[],0)
% sgtitle('Node 10')



%%

% return

figure()
pd1 = histfit(SumIntegExcOut_Vector,[],'beta');
pd1(1).FaceColor = 'none';
pd1(1).EdgeColor = 'none';
pd1(2).Color     = [.2 .2 .2];
xlim([0 0.3])
h2 = xline(0.07,'--r','LineWidth',1.5);
xlabel('Sum Muscle Exc. Cubed')
ylabel('Density')
set(gca,'FontSize',15)
legend([pd1(2),h2],'density','target')
set(gca,'fontsize',16)

% pull out muscle excitations from file
muscle_1_controls = CRBF_excit(time,theta_ref(1:10));
muscle_2_controls = CRBF_excit(time,theta_ref(11:20));
muscle_3_controls = CRBF_excit(time,theta_ref(21:30));
muscle_4_controls = CRBF_excit(time,theta_ref(31:40));
muscle_5_controls = CRBF_excit(time,theta_ref(41:50));
muscle_6_controls = CRBF_excit(time,theta_ref(51:60));


chain_input = vertcat(Draw_Results(:,:,1),Draw_Results(:,:,2),Draw_Results(:,:,3),Draw_Results(:,:,4),Draw_Results(:,:,5));

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


%% main figure

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
xlabel('time (s)')
text(-0.3,1.1,'A','fontsize',10,'fontweight','bold','units','normalized')
ylim([0 100])
yticks([0 50 100])
xlim([0 0.5])
xticks([0 0.25 0.50])
legend([Li,h2],"MCMC","Reference","Location","Northwest",'orientation','horizontal');
legend('boxoff')
box off



subplot(3,3,2)

[Pa,Li] = JackKnife(time_int,rad2deg(velocities_mean),rad2deg(velocities_std),'k',[0.75, 0.75, 0.75]);
set(gca,'fontsize',10)
h2 = plot(time_int,rad2deg(velocity),'r--','LineWidth',1.5);
% h2 = plot(time_int,velocity,'r--','LineWidth',1.5);
ylabel('velocity (deg/s)')
xlabel('time (s)')
text(-0.3,1.1,'B','fontsize',10,'fontweight','bold','units','normalized')
ylim([-600 600])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-600 0 600])
% legend([h1,h2],"Random Draw","Reference","Location","Northwest");
box off

% This section is the muscle excitations cubed
% first set up the prior distribution
mu  = 0.07; % center
sig = 0.05; % width
ma = 0.0;
mi = 0.3;
xp = linspace(max(mi,mu-3*sig),min([ma,mu+3*sig]));
yp = norpf(xp,mu,sig^2);
yn = nordf((mi-mu)/sig)+1-nordf((ma-mu)/sig); % area outside bounds

subplot(3,3,3)

h1 = plot(xp,yp,'--k','LineWidth',1.5);
hold on
pd1 = histfit(SumIntegExcOut_Vector,[],'normal');
pd1(1).FaceColor = 'none';
pd1(1).EdgeColor = 'none';
pd1(2).Color     = [.2 .2 .2];
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

y_max = max(pd1(2).YData);
B = y_max == pd1(2).YData;
index = find([B] == 1);
posterior_center = pd1(2).XData(index); 

integrate = zeros(100,1);
integrate_sum = zeros(100,1);

for j = index:100
        integrate(j-index+1) = (pd1(2).XData(j)-pd1(2).XData(j-1)) * pd1(2).YData(j);
end



% the bottom two rows are for the muscle force trajectories. 
subplot(3,3,4)
[Pa,Li] = JackKnife(time_int,Muscle_1_traj_mean,Muscle_1_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(ReferenceForces(:,1),ReferenceForces(:,2),'r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,1)),'b','LineWidth',2)
xlabel('time (s)')
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
plot(ReferenceForces(:,1),ReferenceForces(:,3),'r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,2)),'b','LineWidth',2)
xlabel('time (s)')
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
h1 = plot(ReferenceForces(:,1),ReferenceForces(:,4),'r','LineWidth',2);
% hold on 
% h2 = plot(Time,(CMCForces_Intp(:,3)),'b','LineWidth',2);
xlabel('time (s)')
ylim([0 625])
yticks([0 300 600])
text(-0.3,1.1,'F','fontsize',10,'fontweight','bold','units','normalized')
box off
set(gca,'fontsize',10)
title('Triceps Med.')
xlim([0 0.5])

subplot(3,3,7)
[Pa,Li] = JackKnife(time_int,Muscle_4_traj_mean,Muscle_4_traj_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(ReferenceForces(:,1),ReferenceForces(:,5),'r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,4)),'b','LineWidth',2)
ylim([0 625])
yticks([0 300 600])
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
plot(ReferenceForces(:,1),ReferenceForces(:,6),'r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,5)),'b','LineWidth',2)
ylim([0 500])
yticks([0 250 500])
legend([Li,h1],'MCMC Mean','Static Optimization','orientation','vertical')
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
plot(ReferenceForces(:,1),ReferenceForces(:,7),'r','LineWidth',2)
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


