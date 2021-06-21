% Run a MCMC to find the plausible muscle forces for a reference motion in
% an elbow musculoskeletal model. This code uses Compact Radial Basis
% functions to create muscle excitation signals to use with the forward
% simulations for each of the six muscles in the model

% written by Russell T Johnson, University of Southern California
% rtjohnso@usc.edu
% Last edited: 5/19/2021

% This depends on following MCMC package for MATLAB: 
% https://mjlaine.github.io/mcmcstat/

% This also depends on having OpenSim 3.3 installed and setting up scripting
% with Matlab :
% https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab

% This code doesn't work on OpenSim 4.x versions directly, there are a
% couple minor coding changes you would have to make, and a few other
% logistical problems that were run into causing the code to take longer
% per iteration. 

% This also uses parallel processing toolbox in Matlab 

clear
clc
close all

user = memory;
memBefore = user.MemUsedMATLAB;

rng(99)

n_pools = 5; 

%% this block loads data from sto file
final_time = 0.5; 
time_interval2 = final_time/100;
time_int = (0:time_interval2:final_time)';

Input_raw = dlmread('arm16_pert4_states.sto','\t',7,0);

position = interp1(Input_raw(:,1),Input_raw(:,2),time_int);
velocity = interp1(Input_raw(:,1),Input_raw(:,3),time_int);

%% Set up Data
data.xdata = time_int;
data.ydata = [time_int,position,velocity];

%% Set up model 
model.ssfun = @Arm16_SS;
model.priorfun = @prior_act;
model.sigma2 = 1;

%% Set up the parameters with "initial values" 

num_parameters = 60;
set_thetamu = 0;
set_thetasig = ones(1,60).*0.2;

Init_Parameter = zeros(num_parameters,n_pools);
Init_Parameter(:,1) = ones(1,60).*-1;

for i = 1:n_pools
    Init_Parameter(:,i) = -5 + 5*rand(1,60);  
end

for i = 1:n_pools
    params(:,i) = {
    % Muscle 1
       {'M10', Init_Parameter(1,i), -30, 30, set_thetamu, set_thetasig(1)} 
       {'M11', Init_Parameter(2,i), -30, 30, set_thetamu, set_thetasig(2)} 
       {'M12', Init_Parameter(3,i), -30, 30, set_thetamu, set_thetasig(3)} 
       {'M13', Init_Parameter(4,i), -30, 30, set_thetamu, set_thetasig(4)} 
       {'M14', Init_Parameter(5,i), -30, 30, set_thetamu, set_thetasig(5)} 
       {'M15', Init_Parameter(6,i), -30, 30, set_thetamu, set_thetasig(6)} 
       {'M16', Init_Parameter(7,i), -30, 30, set_thetamu, set_thetasig(7)} 
       {'M17', Init_Parameter(8,i), -30, 30, set_thetamu, set_thetasig(8)} 
       {'M18', Init_Parameter(9,i), -30, 30, set_thetamu, set_thetasig(9)} 
       {'M19', Init_Parameter(10,i),-30, 30, set_thetamu, set_thetasig(10)} 
    % Muscle 2
       {'M20', Init_Parameter(11,i), -30, 30, set_thetamu, set_thetasig(11)}
       {'M21', Init_Parameter(12,i), -30, 30, set_thetamu, set_thetasig(12)} 
       {'M22', Init_Parameter(13,i), -30, 30, set_thetamu, set_thetasig(13)} 
       {'M23', Init_Parameter(14,i), -30, 30, set_thetamu, set_thetasig(14)} 
       {'M24', Init_Parameter(15,i), -30, 30, set_thetamu, set_thetasig(15)} 
       {'M25', Init_Parameter(16,i), -30, 30, set_thetamu, set_thetasig(16)} 
       {'M26', Init_Parameter(17,i), -30, 30, set_thetamu, set_thetasig(17)} 
       {'M27', Init_Parameter(18,i), -30, 30, set_thetamu, set_thetasig(18)} 
       {'M28', Init_Parameter(19,i), -30, 30, set_thetamu, set_thetasig(19)} 
       {'M29', Init_Parameter(20,i), -30, 30, set_thetamu, set_thetasig(20)} 
    % Muscle 3
       {'M30', Init_Parameter(21,i), -30, 30, set_thetamu, set_thetasig(21)} 
       {'M31', Init_Parameter(22,i), -30, 30, set_thetamu, set_thetasig(22)} 
       {'M32', Init_Parameter(23,i), -30, 30, set_thetamu, set_thetasig(23)} 
       {'M33', Init_Parameter(24,i), -30, 30, set_thetamu, set_thetasig(24)} 
       {'M34', Init_Parameter(25,i), -30, 30, set_thetamu, set_thetasig(25)} 
       {'M35', Init_Parameter(26,i), -30, 30, set_thetamu, set_thetasig(26)} 
       {'M36', Init_Parameter(27,i), -30, 30, set_thetamu, set_thetasig(27)} 
       {'M37', Init_Parameter(28,i), -30, 30, set_thetamu, set_thetasig(28)} 
       {'M38', Init_Parameter(29,i), -30, 30, set_thetamu, set_thetasig(29)} 
       {'M39', Init_Parameter(30,i), -30, 30, set_thetamu, set_thetasig(30)} 
    % Muscle 4
       {'M40', Init_Parameter(31,i), -30, 30, set_thetamu, set_thetasig(31)} % strip off the prior settings and rewrite the prior fun explicitly 
       {'M41', Init_Parameter(32,i), -30, 30, set_thetamu, set_thetasig(32)} 
       {'M42', Init_Parameter(33,i), -30, 30, set_thetamu, set_thetasig(33)} 
       {'M43', Init_Parameter(34,i), -30, 30, set_thetamu, set_thetasig(34)} 
       {'M44', Init_Parameter(35,i), -30, 30, set_thetamu, set_thetasig(35)} 
       {'M45', Init_Parameter(36,i), -30, 30, set_thetamu, set_thetasig(36)} 
       {'M46', Init_Parameter(37,i), -30, 30, set_thetamu, set_thetasig(37)} 
       {'M47', Init_Parameter(38,i), -30, 30, set_thetamu, set_thetasig(38)} 
       {'M48', Init_Parameter(39,i), -30, 30, set_thetamu, set_thetasig(39)} 
       {'M49', Init_Parameter(40,i), -30, 30, set_thetamu, set_thetasig(40)} 
    % Muscle 5
       {'M50', Init_Parameter(41,i), -30, 30, set_thetamu, set_thetasig(41)}
       {'M51', Init_Parameter(42,i), -30, 30, set_thetamu, set_thetasig(42)} 
       {'M52', Init_Parameter(43,i), -30, 30, set_thetamu, set_thetasig(43)} 
       {'M53', Init_Parameter(44,i), -30, 30, set_thetamu, set_thetasig(44)} 
       {'M54', Init_Parameter(45,i), -30, 30, set_thetamu, set_thetasig(45)} 
       {'M55', Init_Parameter(46,i), -30, 30, set_thetamu, set_thetasig(46)} 
       {'M56', Init_Parameter(47,i), -30, 30, set_thetamu, set_thetasig(47)} 
       {'M57', Init_Parameter(48,i), -30, 30, set_thetamu, set_thetasig(48)} 
       {'M58', Init_Parameter(49,i), -30, 30, set_thetamu, set_thetasig(49)} 
       {'M59', Init_Parameter(50,i), -30, 30, set_thetamu, set_thetasig(50)} 
    % Muscle 6
       {'M60', Init_Parameter(51,i), -30, 30, set_thetamu, set_thetasig(51)} 
       {'M61', Init_Parameter(52,i), -30, 30, set_thetamu, set_thetasig(52)} 
       {'M62', Init_Parameter(53,i), -30, 30, set_thetamu, set_thetasig(53)} 
       {'M63', Init_Parameter(54,i), -30, 30, set_thetamu, set_thetasig(54)} 
       {'M64', Init_Parameter(55,i), -30, 30, set_thetamu, set_thetasig(55)} 
       {'M65', Init_Parameter(56,i), -30, 30, set_thetamu, set_thetasig(56)} 
       {'M66', Init_Parameter(57,i), -30, 30, set_thetamu, set_thetasig(57)} 
       {'M67', Init_Parameter(58,i), -30, 30, set_thetamu, set_thetasig(58)} 
       {'M68', Init_Parameter(59,i), -30, 30, set_thetamu, set_thetasig(59)} 
       {'M69', Init_Parameter(60,i), -30, 30, set_thetamu, set_thetasig(60)} 
       
       };
end

%% Run MCMC

% start the clock
tic 

% name of model
nameofModel = 'arm16_millard_rigidtendon.osim';

% settings for MCMC
options.nsimu = 10; % number of iterations
options.waitbar = 0; % don't use the waitbar (doesn't work with parallel anyway)
% define burn-in as a percentage of the number of simulations... 
burn_in = options.nsimu *0.50; 

% must have parallel tool box installed in matlab - check first!
check = contains(struct2array(ver), 'Parallel Computing Toolbox');

if check ~=1
    error('check to see parallel toolbox is installed in matlab')
end

% Open the parallel pools 
poolobj = parpool(n_pools);

% Run the MCMC - Sit Back, this will take a while 
parfor k = 1:n_pools
    [results(:,:,k), chain(:,:,k), s2chain(:,:,k), sschain(:,:,k)]= mcmcrun(model,data,params(:,k),options);
end

% stop the clock
runtime = toc;

% stop the parallel pools
delete(poolobj)

%% Analyze and plot the results

% posterior distribution for each parameter, for each chain
for i = 1:n_pools 
    figure(10+i); clf
    mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel',2);
end

% chain over iterations 
for j = 1:n_pools
    figure(20+j);
    for i =1:size(chain,2)
        subplot(6,10,i)
        plot(chain(:,i,j))
    %     ylim([ 0 1]);
        set(gca,'fontsize',14)
        if i == 1
            ylabel('amplitude')
            title('Tri long')
        end
        if i == 11
            ylabel('amplitude')
            title('Tri lat')
        end
        if i == 21
            ylabel('amplitude')
            title('Tri med')
        end
        if i == 31
            ylabel('amplitude')
            title('Biceps LH')
        end
        if i == 41
            ylabel('amplitude')
            title('Biceps SH')
        end
        if i == 51
            ylabel('amplitude')
            title('Brachialis')
        end

    end
end

% all the chains for all 5 chains - would need to manually change with
% different number of parallel pools 
figure(29)
for i =1:size(chain,2)
    subplot(6,10,i)
    plot(chain(:,i,1),'r','LineWidth',2)
    hold on 
    plot(chain(:,i,2),'b','LineWidth',2)
    hold on 
    plot(chain(:,i,3),'g','LineWidth',2)
    hold on 
    plot(chain(:,i,4),'k','LineWidth',2)
    hold on 
    plot(chain(:,i,5),'c','LineWidth',2)
    xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.2);
end

% Create the rank plots to understand convergence 
rank_plot_Arm16_10CRBVs_fxn(chain(burn_in:end,:,:))

%% Choose 1 chain and then plot some information from that one
% Random draw for a chain to evaluate performance

choose = 1; % chose chain
chain_plot = chain(:,:,choose);

for k = 1:20
   draw(k) = randi([burn_in options.nsimu]);
   Draw_Results(k,:) = chain_plot(draw(k),:);
   yy = Arm16_SimManager_controls_CRBF_6musc(Draw_Results(k,:));
   positions_draw(:,k) = interp1(yy(:,1),yy(:,2),time_int);
   velocities_draw(:,k) = interp1(yy(:,1),yy(:,3),time_int);
end

for i = 1:10 % this 10 is for the number of nodes per muscle. not the number of random draws.. 
    Muscle_draw1(:,i) = Draw_Results(:,i);
    Muscle_draw2(:,i) = Draw_Results(:,i+10);
    Muscle_draw3(:,i) = Draw_Results(:,i+20);
    Muscle_draw4(:,i) = Draw_Results(:,i+30);
    Muscle_draw5(:,i) = Draw_Results(:,i+40);
    Muscle_draw6(:,i) = Draw_Results(:,i+50);
end

tfinal = 0.5;
time = (0:0.001:tfinal)';
for i = 1:20
    Muscle_draw1_traj(:,i) = CRBF_excit(time,Muscle_draw1(i,:));
    Muscle_draw2_traj(:,i) = CRBF_excit(time,Muscle_draw2(i,:));
    Muscle_draw3_traj(:,i) = CRBF_excit(time,Muscle_draw3(i,:));    
    Muscle_draw4_traj(:,i) = CRBF_excit(time,Muscle_draw4(i,:));
    Muscle_draw5_traj(:,i) = CRBF_excit(time,Muscle_draw5(i,:));
    Muscle_draw6_traj(:,i) = CRBF_excit(time,Muscle_draw6(i,:));
end

for i = 1:20
    controls(:,1,i) = CRBF_excit(time,Muscle_draw1(i,:));
    controls(:,2,i) = CRBF_excit(time,Muscle_draw2(i,:));
    controls(:,3,i) = CRBF_excit(time,Muscle_draw3(i,:));    
    controls(:,4,i) = CRBF_excit(time,Muscle_draw4(i,:));
    controls(:,5,i) = CRBF_excit(time,Muscle_draw5(i,:));
    controls(:,6,i) = CRBF_excit(time,Muscle_draw6(i,:));
end
SumIntegExcOut = zeros(20,1);
h = 0.001;

for i = 1:20 
    for k = 1:6
        SumIntegExcOut(i) = SumIntegExcOut(i) + (h * trapz(controls(:,k,i).^3));
    end
end

% muscle excitations for each of the six muscles from the random draws 
figure(52)
subplot(2,3,1)
for i = 1:20
    plot(time,Muscle_draw1_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('tri long')
ylabel('excitation')
ylim([0 1])

subplot(2,3,2)
for i = 1:20
    plot(time,Muscle_draw2_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('tri lat')
ylabel('excitation')
ylim([0 1])

subplot(2,3,3)
for i = 1:20
    plot(time,Muscle_draw3_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('tri med')
xlabel('time (s)')
ylabel('excitation')
ylim([0 1])


subplot(2,3,4)
for i = 1:20
    plot(time,Muscle_draw4_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('biceps lh')
ylabel('excitation')
ylim([0 1])

subplot(2,3,5)
for i = 1:20
    plot(time,Muscle_draw5_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('biceps sh')
ylabel('excitation')
ylim([0 1])

subplot(2,3,6)
for i = 1:20
    plot(time,Muscle_draw6_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('brachialis')
xlabel('time (s)')
ylabel('excitation')
ylim([0 1])

% the sum of squared errors from the results, tracking how well the motion
% matches the reference data 
figure(53)
plot(sschain(:,:,choose),'color', [0.4660 0.6740 0.1880])
ylabel('ss value')
xlabel('iteration')

% plot the kinematics (positions and velocities) from the random draws
% versus the reference data 

figure(81)
subplot(2,1,1)
for j = 1:10
    h1 = plot(time_int,positions_draw(:,j),'color',[.17 .17 .17],'LineWidth',2);
    hold on 
end
set(gca,'fontsize',16)
h2 = plot(time_int,position,'r--','LineWidth',1.5);
ylabel('position (rad)')
xlabel('time (s)')
legend([h1,h2],"Random Draw","Reference","Location","Northwest");
box off

subplot(2,1,2)
for j = 1:10
    plot(time_int,velocities_draw(:,j),'color',[.17 .17 .17],'LineWidth',2)
    hold on 
end
set(gca,'fontsize',16)
plot(time_int,velocity,'r--','LineWidth',1.5)
xlabel('time (s)')
ylabel('velocity (rad/s)')

% Plot the correlations within a muscle - uncomment if you want to see the
% correlation plots 
% names = results_one.names;

% figure(101)
% pairs_v3(chain(burn_in:end,1:10,choose),'panellims',names(1:10),[],0)
% title('tri long')
% 
% figure(102)
% pairs_v3(chain(burn_in:end,11:20,choose),'panellims',names(11:20),[],0)
% title('tri lat')
% 
% figure(103)
% pairs_v3(chain(burn_in:end,21:30,choose),'panellims',names(21:30),[],0)
% title('tri med')
% 
% figure(104)
% pairs_v3(chain(burn_in:end,31:40,choose),'panellims',names(31:40),[],0)
% title('biceps lh')
% 
% figure(105)
% pairs_v3(chain(burn_in:end,41:50,choose),'panellims',names(41:50),[],0)
% title('biceps sh')
% 
% figure(106)
% pairs_v3(chain(burn_in:end,51:60,choose),'panellims',names(51:60),[],0)
% title('brachiorad')
% 
% %  Plot the correlations within a node
% figure(111)
% pairs_v3(chain(burn_in:end,1:10:60,choose),'panellims',names(1:10:60),[],0)
% sgtitle('Node 1')
% 
% figure(112)
% pairs_v3(chain(burn_in:end,2:10:60,choose),'panellims',names(2:10:60),[],0)
% title('Node 2')
% 
% figure(113)
% pairs_v3(chain(burn_in:end,3:10:60,choose),'panellims',names(3:10:60),[],0)
% title('Node 3')
% 
% figure(114)
% pairs_v3(chain(burn_in:end,4:10:60,choose),'panellims',names(4:10:60),[],0)
% title('Node 4')
% 
% figure(115)
% pairs_v3(chain(burn_in:end,5:10:60,choose),'panellims',names(5:10:60),[],0)
% title('Node 5')
% 
% figure(116)
% pairs_v3(chain(burn_in:end,6:10:60,choose),'panellims',names(6:10:60),[],0)
% title('Node 6')
% 
% figure(117)
% pairs_v3(chain(burn_in:end,7:10:60,choose),'panellims',names(7:10:60),[],0)
% title('Node 7')
% 
% figure(118)
% pairs_v3(chain(burn_in:end,8:10:60,choose),'panellims',names(8:10:60),[],0)
% title('Node 8')
% 
% figure(119)
% pairs_v3(chain(burn_in:end,9:10:60,choose),'panellims',names(9:10:60),[],0)
% title('Node 9')
% 
% figure(120)
% pairs_v3(chain(burn_in:end,10:10:60,choose),'panellims',names(10:10:60),[],0)
% title('Node 10')

filename = ['chain_results_',datestr(now,'yyyymmddTHHMMSS'),'.mat'];

save(filename,'chain');

user = memory;
memAfter = user.MemUsedMATLAB;

% display some stuff about the MCMC run
disp(['========================== ' ])
disp(['elapsed time = ' num2str(runtime/60) ' min'])
disp(['change in Matlab memory use = ' num2str((memAfter - memBefore)/1e6) ' MB'])
disp('   ')

%% End of main script - below are functions for the posterior probability and prior

function ss = Arm16_SS(theta,data)

    ydata_pos = data.ydata(:,2); % this is for only caring about position
    ydata_vel = data.ydata(:,3);
       
    y = Arm16_fun(theta); % run the forward simulation
    
    ymodel_pos = y(:,1); % position
    ymodel_vel = y(:,2); % velocity
   
    likelihood(1) = sum(((ymodel_pos - ydata_pos).^2)/0.05^2);
    likelihood(2) = sum(((ymodel_vel - ydata_vel).^2)/0.50^2);

    ss = sum(likelihood); 
end

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
    
    P = 0.07; % sum of muscle excitations from reference motion
    prior = ((SumIntegExc - P)/.05)^2;
end 

function y=Arm16_fun(theta)
    % code to run forward simulation
    time_100 = [0:0.005:0.5];
    raw_out = Arm16_SimManager_controls_CRBF_6musc(theta); 
    y(:,1) =  interp1(raw_out(:,1),raw_out(:,2),time_100);
    y(:,2) =  interp1(raw_out(:,1),raw_out(:,3),time_100);

end