clear
clc
close all

user = memory;
memBefore = user.MemUsedMATLAB;

% global trial 
% 
% trial = 1; 

% run a mcmc algorithm to predict/recover the mass, stiffness, and damping values
% for a mass-spring-damper system based on some data

% The update to this code (version_2, Chebys) uses Compact Radial Basis Fxns to
% create muscle excitation signals to use with the forward simulations. 

% This next update incorporates all six muscles of the arm 16 model.
% Including the antagonists muscles will help with co-contraction measures

% Import the OpenSim modeling classes
% import org.opensim.modeling.*

% Declare osimModel osimState as global
% global osimModel osimState;
% global trial 
% trial = trial+1;

% This is for opensim 4.1 but 4.1 is still dumb, need to figure out how to
% supress output to command window
% path='C:\OpenSim 4.1\Geometry';
% ModelVisualizer.addDirToGeometrySearchPaths(path);

% Path name of OsimModel -- !! When you change it here, you need to change
% it on the *controls*.m script too !!!!
n_pools = 7; 


% Read in the osim model
% osimModel = Model(nameofModel);

%% this block loads data from sto file
final_time = 0.5; 
time_interval2 = final_time/100;
time_int = (0:time_interval2:final_time)';

% Input_raw = dlmread('arm16_pert4_states.sto','\t',7,0);
input_struct = load('arm16_noise_kinematics.mat');
input_controls_ref = readmatrix('arm16_Tracking_p50_degroote_w75_cubed_v6controls.sto','FileType','text');
input_data = input_struct.new_noise_data; 

position = interp1(input_data(:,1),input_data(:,2),time_int);
velocity = interp1(input_data(:,1),input_data(:,3),time_int);


%% Set up Data
data.xdata = time_int;
data.ydata = [time_int,position,velocity];

%% Set up model 
% return

model.ssfun = @Arm16_SS;
% model.priorfun = @(th,alpha,beta) sum(1./(betapdf(th,alpha,beta))); %
% This is my attempt to use a beta distribution for a prior with shape
% factors A=1 and B= 3; model.priorfun = @(th) Cheby_excit_prior(th); %
% model.priorfun = @Cheby_excit_prior; %

% HACK the system and make it so that the prior always returns ZERO.
% Instead we are just going to make the liklihood all encapsulating...
% model.priorfun = @zero_prior;
% model.priorfun = @prior_act;
% uniform distribution for a prior with min and max values - I'm not sure
% it needs to be in the demoninator or not - may be academic... 
% model.priortype = 2; 
model.sigma2 = 1;

%% Set up the parameters with "initial values" 

num_parameters = 60;
set_thetamu = -15;
set_thetasig = ones(1,60).*5;

% For beta prior, shape parameters are alpha and beta
% alpha = 1;
% beta = 3; 

% Create initial parameters from beta distribution, alpha =1; beta =3;
% Init_Parameter = betarnd(1,3,num_parameters,1);

rng(111)

Init_Parameter = zeros(num_parameters,n_pools);

% for i = 1:n_pools
%     Init_Parameter(1:30,i) = (ones(30,1).*-1.6)+(randn(30,1).*0.03);
% % These inital parameters are from a CMC fit!
% %     Init_Parameter(25:32,i) = [-2.66, -0.20, -2.64, -1.62, -0.44, -0.95, -1.0, -0.04];
% %     Init_Parameter(33:40,i) = [-1.62, -1.52, -1.27, -0.44, -0.87, -1.45, -0.64, -0.67];
% %     Init_Parameter(41:48,i) = [ -1.19, -1.16, 0.18, -1.25, -0.66, -1.02, -0.92, -0.65];
%     Init_Parameter(31:60,i) = (ones(30,1).*-1.0)+(randn(30,1).*0.03);
%     
% %     Init_Parameter(1:24,i) = (ones(24,1).*-2)+(randn(24,1).*0.05);
% %     Init_Parameter(25:48,i) = (ones(24,1).*-0.85)+(randn(24,1).*0.05);
% end
% %   {'Name', Initial Proposal, min, max, prior_mu,   prior_sig}


% Inits(1:10)  = [-4.14420, -0.402106, -3.97480, -0.51432, -2.8496,  0.38420, -2.2335, -2.8773, -1.45525, -2.78680];
% Inits(11:20) = [-4.52530, -0.085629, -4.31510, -0.18319, -2.1970,  0.97354, -2.1090, -3.2583, -1.06420, -3.37399];
% Inits(21:30) = [-4.47320, -0.159957, -4.16200, -0.25489, -2.1679,  0.90320, -2.1217, -3.1114, -1.20240, -3.21812];
% Inits(31:40) = [ 0.71920, -5.424749,  2.59178, -5.21709, -1.6753, -2.37993, -2.3760, -0.4516, -3.35891, -0.94957];
% Inits(41:50) = [-0.39960, -2.692176,  1.50546, -3.20474, -90.63, -0.14599, -4.1915,  1.3608, -3.79824, -0.70829]; 
% Inits(51:60) = [-1.38806,  0.035885,  0.60393, -2.34746, -90.06, -2.95874, -1.6006, -1.4812,  0.73418, -3.27991]; 
% 
% for i = 1: n_pools
%     Init_Parameter(:,i) = Inits + rand(num_parameters,1).*0.05; 
% end

% Inits = theta_ref; 
% 
% for i = 1:length(Inits) % Updating the bounds on min/max, so checking that the OG points don't start outside those bounds
%     if Inits(i) < -25
%         Inits(i) = -25+randn(1,1);
%     end    
% end
% Scale_Factor = 0.03;
% set_thetasig = abs(Inits.*0.015);
% for i = 1:length(set_thetasig)
%     if set_thetasig(i) < 0.01
%         set_thetasig(i) = 0.01;
%     end
% end

Init_Parameter(:,1) = ones(1,60).*-1;


% load('matrix_crbf_amps.mat')
% Init_Parameter(:,1) = vertcat(Matrix(:,1),Matrix(:,2),Matrix(:,3),Matrix(:,4),Matrix(:,5),Matrix(:,6));

for i = 1:n_pools
    Init_Parameter(:,i) = -15 + 10*rand(1,60);  
end

Init_Parameter(61,:) = 0; 

% return

for i = 1:n_pools
    params(:,i) = {
    % Muscle 1
       {'M10', Init_Parameter(1,i), -30, 30, set_thetamu, set_thetasig(1)} % strip off the prior settings and rewrite the prior fun explicitly 
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

       
       {'pretemp', Init_Parameter(61,i), -Inf, Inf, 0, 0.5}
       
%        {'p_init',0,-1,1,0,0.05}
%        {'v_init',0,-5,5,0,0.10}
       
       };
end




%% First generate an initial chain.
% rng(42);
tic 

% global osimModel osimState;
% poolobj = parpool(n_pools);

% nameofModel = Model('arm16_millard_rigidtendon_tris.osim');
nameofModel = 'arm16_millard_rigidtendon.osim';

% parfor i = 1:n_pools % start parallel for-loop and "assign" osimModel & osimState for each worker
%      startParallelCoreComputing(nameofModel);
% end

options.nsimu = 25000;
% options.stats = 1; 
% options.stats2 = 1; 
options.waitbar = 0;

% load old results 
% load results_20220103T095926
% old_results = results; 
% clear results 
% trying these to see what happens
% options.burnintime = options.nsimu * 0.2; 

% Open the parallel pools 
poolobj = parpool(n_pools);

parfor k = 1:n_pools
    [results(:,:,k), chain(:,:,k), s2chain(:,:,k), sschain(:,:,k)]= mcmcrun(model,data,params(:,k),options);
end

% if restarting from old results 
% parfor k = 1:n_pools
%     [results(:,:,k), chain(:,:,k), s2chain(:,:,k), sschain(:,:,k)]= mcmcrun(model,data,params(:,k),options,old_results(:,:,k));
% end

runtime = toc;

% stop the parallel pools
delete(poolobj)

% figure(93); clf
% mcmcplot(chain,[],results,'denspanel',2);


%% Analyze and plot the results


% define burn-in as a percentage of the number of simulations... 
burn_in = options.nsimu *0.50;

% results_one = results(:,:,1);
% results_two = results(:,:,2);
% results_thr = results(:,:,3); 
% results_for = results(:,:,4);

% return

for i = 1:n_pools 
    figure(10+i); clf
    mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel',2);
end

for j = 1:n_pools
    figure(20+j);
    for i =1:60
        subplot(6,10,i)
        plot(chain(:,i,j))
    %     ylim([ 0 1]);
        set(gca,'fontsize',14)
        if i == 1
            ylabel('excitation')
            title('Tri long')
        end
        if i == 11
            ylabel('excitation')
            title('Tri lat')
        end
        if i == 21
            ylabel('excitation')
            title('Tri med')
        end
        if i == 31
            ylabel('excitation')
            title('Biceps LH')
        end
        if i == 41
            ylabel('excitation')
            title('Biceps SH')
        end
        if i == 51
            ylabel('excitation')
            title('Brachialis')
        end

    end
end

% plot temperature parameter 

      offset = 4; % try between 2 and 4
      amplitude = 0.4; % try between 0 and 1
      
for k = 1:n_pools
    for j = 1:size(chain,1)
        denominator(j,k) = temperature(chain(j,61,k),offset,amplitude);
    end
end

for k = 1:n_pools
    for j = 1:size(chain,1)
        prior2(j,k) = (1/2)*(chain(j,61,k))^2;
    end
end

figure(299)
subplot(3,1,1)
plot(chain(:,61,1),'color','#0072BD','LineWidth',2)
hold on 
plot(chain(:,61,2),'color','#D95319','LineWidth',2)
hold on 
plot(chain(:,61,3),'color','#EDB120','LineWidth',2)
hold on 
plot(chain(:,61,4),'color','#7E2F8E','LineWidth',2)
hold on 
plot(chain(:,61,5),'color','#77AC30','LineWidth',2)
hold on 
plot(chain(:,61,6),'color','#4DBEEE','LineWidth',2)
hold on 
plot(chain(:,61,7),'color','#A2142F','LineWidth',2)   
xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.2);
title('pre temp')

subplot(3,1,2)
plot(denominator(:,1),'color','#0072BD','LineWidth',2)
hold on 
plot(denominator(:,2),'color','#D95319','LineWidth',2)
hold on 
plot(denominator(:,3),'color','#EDB120','LineWidth',2)
hold on 
plot(denominator(:,4),'color','#7E2F8E','LineWidth',2)
hold on 
plot(denominator(:,5),'color','#77AC30','LineWidth',2)
hold on 
plot(denominator(:,6),'color','#4DBEEE','LineWidth',2)
hold on 
plot(denominator(:,7),'color','#A2142F','LineWidth',2)   
xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.2);
title('temperature')

subplot(3,1,3)
plot(prior2(:,1),'color','#0072BD','LineWidth',2)
hold on 
plot(prior2(:,2),'color','#D95319','LineWidth',2)
hold on 
plot(prior2(:,3),'color','#EDB120','LineWidth',2)
hold on 
plot(prior2(:,4),'color','#7E2F8E','LineWidth',2)
hold on 
plot(prior2(:,5),'color','#77AC30','LineWidth',2)
hold on 
plot(prior2(:,6),'color','#4DBEEE','LineWidth',2)
hold on 
plot(prior2(:,7),'color','#A2142F','LineWidth',2)   
xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.2);
title('prior2')
%     
figure(33)
for i =1:60
    subplot(6,10,i)
    plot(chain(:,i,1),'color','#0072BD','LineWidth',2)
    hold on 
    plot(chain(:,i,2),'color','#D95319','LineWidth',2)
    hold on 
    plot(chain(:,i,3),'color','#EDB120','LineWidth',2)
    hold on 
    plot(chain(:,i,4),'color','#7E2F8E','LineWidth',2)
    hold on 
    plot(chain(:,i,5),'color','#77AC30','LineWidth',2)
    hold on 
    plot(chain(:,i,6),'color','#4DBEEE','LineWidth',2)
    hold on 
    plot(chain(:,i,7),'color','#A2142F','LineWidth',2)   
    xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.2);
end
 legend('c1','c2','c3','c4','c5','c6','c7','orientation','horizontal')
 legend('boxoff')
 
% try pushing the chain through the log it transform equation and
% recalculate the rank plots? 
% for i = 1:size(chain,1)
%     for j = 1:size(chain,2)
%         for k = 1:size(chain,3) 
%             if chain(i,j,k)>0 
%                 chain_transform(i,j,k) = 1./(1+exp(-chain(i,j,k)));
%             else
%                 chain_transform(i,j,k) = exp(chain(i,j,k))./(1+exp(chain(i,j,k)));
%             end
%         end
%     end
% end


% Create the rank plots to check for convergence 
rank_plot_Arm16_10CRBVs_fxn(chain(burn_in:end,:,:))

% Output = chainstats(chain,results);
% 
% Means = Output(:,1);
% Stds  = Output(:,2);
%% Choose 1 chain and then plot that one 
% Random draw for a chain to evaluate performance

choose = 2; % chose chain/result 1

chain_plot = chain(:,:,choose);
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

color_scheme = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};

figure()

set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 8])

histogram(SumIntegExcOut(:,1),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{1},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,2),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{2},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,3),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{3},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,4),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{4},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,5),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{5},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,6),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{6},'FaceAlpha',0.5);
hold on 
histogram(SumIntegExcOut(:,7),'NumBins',40,'BinLimits',[0 1],'FaceColor',color_scheme{7},'FaceAlpha',0.5);
xlabel('Sum Muscle Excitation Cubed')
set(gca,'FontSize',10)
% ylim([0 8])
% yticks([0 4 8])
box off
% xticks([0 0.1 0.2 0.3])
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','location','northeast')
legend('boxoff')

SumIntegExcOut_Vector = vertcat(SumIntegExcOut(:,1),SumIntegExcOut(:,2),SumIntegExcOut(:,3),SumIntegExcOut(:,4),SumIntegExcOut(:,5),SumIntegExcOut(:,6),SumIntegExcOut(:,7));


chain_input = vertcat(Draw_Results(:,:,1),Draw_Results(:,:,2),Draw_Results(:,:,3),Draw_Results(:,:,4),Draw_Results(:,:,5),Draw_Results(:,:,6),Draw_Results(:,:,7));

for i = 1:size(chain_input,1)
    [kinematics_out(:,:,i),force_out(:,:,i),controls_out(:,:,i)] = Arm16_SimManager_controls_CRBF_6musc_wForce(chain_input(i,:));
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

Muscle_1_cont_mean = mean(controls_out(:,1,:),3);
Muscle_2_cont_mean = mean(controls_out(:,2,:),3);
Muscle_3_cont_mean = mean(controls_out(:,3,:),3);
Muscle_4_cont_mean = mean(controls_out(:,4,:),3);
Muscle_5_cont_mean = mean(controls_out(:,5,:),3);
Muscle_6_cont_mean = mean(controls_out(:,6,:),3);

Muscle_1_cont_std(:,1) = std(controls_out(:,1,:),0,[3]);
Muscle_2_cont_std(:,1) = std(controls_out(:,2,:),0,[3]);
Muscle_3_cont_std(:,1) = std(controls_out(:,3,:),0,[3]);
Muscle_4_cont_std(:,1) = std(controls_out(:,4,:),0,[3]);
Muscle_5_cont_std(:,1) = std(controls_out(:,5,:),0,[3]);
Muscle_6_cont_std(:,1) = std(controls_out(:,6,:),0,[3]);




% Create the rank plots to check for convergence 
% rank_plot_Arm16_10CRBVs_fxn(chain(burn_in:end,:,:))
% muscles = {'Tri Long', 'Tri Lat','Tri Med','Biceps LH', 'Biceps SH','Brachior'};
% muscle 4 (Biceps LH)
% nNodes = 10;
% Chain_rank_plot_muscle_10node(chain(:,((1*nNodes)-9):(1*nNodes),:),'Triceps Long',burn_in)
% Chain_rank_plot_muscle_10node(chain(:,((2*nNodes)-9):(2*nNodes),:),'Triceps Lat',burn_in)
% Chain_rank_plot_muscle_10node(chain(:,((3*nNodes)-9):(3*nNodes),:),'Triceps Med',burn_in)
% Chain_rank_plot_muscle_10node(chain(:,((4*nNodes)-9):(4*nNodes),:),'Biceps LH',burn_in)
% Chain_rank_plot_muscle_10node(chain(:,((5*nNodes)-9):(5*nNodes),:),'Biceps SH',burn_in)
% Chain_rank_plot_muscle_10node(chain(:,((6*nNodes)-9):(6*nNodes),:),'Brachialis',burn_in)

figure(61)

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
ylim([-600 600])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-600 0 600])
% legend([h1,h2],"Random Draw","Reference","Location","Northwest");
box off

% This section is the muscle excitations cubed
% first set up the prior distribution
mu  = 0.10; % center
sig = 0.05; % width
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
xlim([0 1])
% ylim([0 40])
xlabel('Sum Muscle Exc. Cubed')
ylabel('density')
text(-0.3,1.1,'C','fontsize',10,'fontweight','bold','units','normalized')
% legend([pd(2),h2],'density','target')
legend('boxoff')
set(gca,'fontsize',10)
box off



% calculate the width of posterior: 

% y_max = max(pd1(2).YData);
% B = y_max == pd1(2).YData;
% index = find([B] == 1);
% posterior_center = pd1(2).XData(index); 
% 
% integrate = zeros(100,1);
% integrate_sum = zeros(100,1);
% 
% for j = index:100
%         integrate(j-index+1) = (pd1(2).XData(j)-pd1(2).XData(j-1)) * pd1(2).YData(j);
% end

% integrate_sum(1) = integrate(1);
% for k = 2:100
%     integrate_sum(k) = integrate_sum(k-1)+integrate(k);
% end
% C = integrate_sum<0.341;
% width = index + sum(C);
% posterior_std = pd(2).XData(width)-posterior_center;
% 
% posterior_center
% posterior_std
% 

% the bottom two rows are for the muscle force trajectories. 
subplot(3,3,4)
[Pa,Li] = JackKnife(time_int,Muscle_1_traj_mean,Muscle_1_traj_std,'k',[0.75, 0.75, 0.75]);
% hold on 
% plot(StatOpForces(:,1),StatOpForces(:,2),'--r','LineWidth',2)
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
% hold on 
% plot(StatOpForces(:,1),StatOpForces(:,3),'--r','LineWidth',2)
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
% hold on 
% h1 = plot(StatOpForces(:,1),StatOpForces(:,4),'--r','LineWidth',2);
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
% hold on 
% plot(StatOpForces(:,1),StatOpForces(:,5),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,4)),'b','LineWidth',2)
ylim([0 625])
yticks([0 300 600])
% legend([Pa],'');
% legend boxoff
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
% hold on 
% plot(StatOpForces(:,1),StatOpForces(:,6),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,5)),'b','LineWidth',2)
ylim([0 500])
yticks([0 250 500])
% legend([Li,h1],'MCMC Mean +/- SD','Reference','orientation','vertical');
% labelhandles(4).FaceColor = [0.75,0.75,0.75];
% labelhandles(3).YData = [0.83 0.83];
% labelhandles(4).XData = [0.0460 0.25]; labelhandles(4).YData = [0.760 0.760];
% legend boxoff
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
% hold on 
% plot(StatOpForces(:,1),StatOpForces(:,7),'--r','LineWidth',2)
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


figure(65)

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
ylim([-600 600])
xlim([0 0.5])
xticks([0 0.25 0.50])
yticks([-600 0 600])
% legend([h1,h2],"Random Draw","Reference","Location","Northwest");
box off

% This section is the muscle excitations cubed
% first set up the prior distribution
mu  = 0.00; % center
sig = 0.10; % width
ma = 0.0;
mi = 1;
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
xlim([0 1])
% ylim([0 40])
xlabel('Sum Muscle Exc. Cubed')
ylabel('density')
text(-0.3,1.1,'C','fontsize',10,'fontweight','bold','units','normalized')
% legend([pd(2),h2],'density','target')
legend('boxoff')
set(gca,'fontsize',10)
box off



% calculate the width of posterior: 

% y_max = max(pd1(2).YData);
% B = y_max == pd1(2).YData;
% index = find([B] == 1);
% posterior_center = pd1(2).XData(index); 
% 
% integrate = zeros(100,1);
% integrate_sum = zeros(100,1);
% 
% for j = index:100
%         integrate(j-index+1) = (pd1(2).XData(j)-pd1(2).XData(j-1)) * pd1(2).YData(j);
% end

% integrate_sum(1) = integrate(1);
% for k = 2:100
%     integrate_sum(k) = integrate_sum(k-1)+integrate(k);
% end
% C = integrate_sum<0.341;
% width = index + sum(C);
% posterior_std = pd(2).XData(width)-posterior_center;
% 
% posterior_center
% posterior_std
% 

% the bottom two rows are for the muscle force trajectories. 
subplot(3,3,4)
[Pa,Li] = JackKnife(time_int,Muscle_1_cont_mean,Muscle_1_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(input_controls_ref(:,1),input_controls_ref(:,2),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,1)),'b','LineWidth',2)
% xlabel('time (s)')
box off
set(gca,'fontsize',10)
title('Triceps Long.')
text(-0.3,1.1,'D','fontsize',10,'fontweight','bold','units','normalized')
ylabel('Excitation')
ylim([0 1])
yticks([0 0.5 1])
xlim([0 0.5])
xticks([0 0.25 0.50])

% -----
subplot(3,3,5)
[Pa,Li] = JackKnife(time_int,Muscle_2_cont_mean,Muscle_2_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(input_controls_ref(:,1),input_controls_ref(:,3),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,2)),'b','LineWidth',2)
% xlabel('time (s)')
box off
set(gca,'fontsize',10)
title('Triceps Lat.')
text(-0.3,1,'E','fontsize',10,'fontweight','bold','units','normalized')
% ylabel('Force (N)')
ylim([0 1])
yticks([0 0.5 1])
xlim([0 0.5])
xticks([0 0.25 0.50])


subplot(3,3,6)
[Pa,Li] = JackKnife(time_int,Muscle_3_cont_mean,Muscle_3_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
h1 = plot(input_controls_ref(:,1),input_controls_ref(:,4),'--r','LineWidth',2);
% hold on 
% h2 = plot(Time,(CMCForces_Intp(:,3)),'b','LineWidth',2);
% xlabel('time (s)')
ylim([0 1])
yticks([0 0.5 1])
text(-0.3,1.1,'F','fontsize',10,'fontweight','bold','units','normalized')
box off
set(gca,'fontsize',10)
title('Triceps Med.')
xlim([0 0.5])
xticks([0 0.25 0.50])

subplot(3,3,7)
[Pa,Li] = JackKnife(time_int,Muscle_4_cont_mean,Muscle_4_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(input_controls_ref(:,1),input_controls_ref(:,5),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,4)),'b','LineWidth',2)
ylim([0 1])
yticks([0 0.5 1])
% legend([Pa],'');
% legend boxoff
box off
set(gca,'fontsize',10)
text(-0.3,1.1,'G','fontsize',10,'fontweight','bold','units','normalized')
title('Biceps LH')
ylabel('Excitation')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])


subplot(3,3,8)
[Pa,Li] = JackKnife(time_int,Muscle_5_cont_mean,Muscle_5_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(input_controls_ref(:,1),input_controls_ref(:,6),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,5)),'b','LineWidth',2)
ylim([0 1])
yticks([0 0.5 1])
% legend([Li,h1],'MCMC Mean +/- SD','Reference','orientation','vertical');
% labelhandles(4).FaceColor = [0.75,0.75,0.75];
% labelhandles(3).YData = [0.83 0.83];
% labelhandles(4).XData = [0.0460 0.25]; labelhandles(4).YData = [0.760 0.760];
% legend boxoff
box off
set(gca,'fontsize',10)
title('Biceps SH')
text(-0.3,1.1,'H','fontsize',10,'fontweight','bold','units','normalized')
% ylabel('Force (N)')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])

subplot(3,3,9)
[Pa,Li] = JackKnife(time_int,Muscle_6_cont_mean,Muscle_6_cont_std,'k',[0.75, 0.75, 0.75]);
hold on 
plot(input_controls_ref(:,1),input_controls_ref(:,7),'--r','LineWidth',2)
% hold on 
% plot(Time,(CMCForces_Intp(:,6)),'b','LineWidth',2)
text(-0.3,1.1,'I','fontsize',10,'fontweight','bold','units','normalized')
ylim([0 1])
yticks([0 0.5 1])
box off
set(gca,'fontsize',10)
title('Brachialis')
% ylabel('Force (N)')
xlabel('time (s)')
xlim([0 0.5])
xticks([0 0.25 0.50])


% time_plot = 0:0.005:0.5;

% pull out muscle excitations from file
% muscle_1_controls = CRBF_excit(time,theta_ref(1:10));
% muscle_2_controls = CRBF_excit(time,theta_ref(11:20));
% muscle_3_controls = CRBF_excit(time,theta_ref(21:30));
% muscle_4_controls = CRBF_excit(time,theta_ref(31:40));
% muscle_5_controls = CRBF_excit(time,theta_ref(41:50));
% muscle_6_controls = CRBF_excit(time,theta_ref(51:60));

% Need to rewrite this to update with chebys
% contourplot_MCMC(chain,time,time_int,controls,[],[],[],[])
% contourplot_MCMC(chain,time,time_int,controls)

figure(52)
subplot(2,3,1)
% [fillhandle,msg]=jbfill(time,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time,muscle_1_controls,'b','LineWidth',2)
% hold on 
for i = 1:20
    plot(time,Muscle_draw1_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('tri long')
ylabel('excitation')
ylim([0 1])

subplot(2,3,2)
% [fillhandle,msg]=jbfill(time,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time,muscle_2_controls,'b','LineWidth',2)
% hold on 
for i = 1:20
    plot(time,Muscle_draw2_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('tri lat')
ylabel('excitation')
ylim([0 1])

subplot(2,3,3)
% [fillhandle,msg]=jbfill(time,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time,muscle_3_controls,'b','LineWidth',2)
% hold on 
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
% [fillhandle,msg]=jbfill(time,Muscle_4_exc_upper',Muscle_4_exc_lower');
% plot(time,muscle_4_controls,'b','LineWidth',2)
% hold on 
for i = 1:20
    plot(time,Muscle_draw4_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('biceps lh')
ylabel('excitation')
ylim([0 1])

subplot(2,3,5)
% [fillhandle,msg]=jbfill(time,Muscle_5_exc_upper',Muscle_5_exc_lower');
% plot(time,muscle_5_controls,'b','LineWidth',2)
% hold on 
for i = 1:20
    plot(time,Muscle_draw5_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('biceps sh')
ylabel('excitation')
ylim([0 1])

subplot(2,3,6)
% [fillhandle,msg]=jbfill(time,Muscle_6_exc_upper',Muscle_6_exc_lower');
% plot(time,muscle_6_controls,'b','LineWidth',2)
% hold on 
for i = 1:20
    plot(time,Muscle_draw6_traj(:,i),'color',[.17 .17 .17],'LineWidth',1.5)
    hold on 
end
set(gca,'fontsize',16)
title('brachialis')
xlabel('time (s)')
ylabel('excitation')
ylim([0 1])



figure(53)
for i = 1:n_pools
    plot(sschain(:,:,i),'color', color_scheme{i},'LineWidth',2)
    hold on 
end
ylabel('ss value')
xlabel('iteration')
legend('Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','Chain 6','Chain 7','location','northeast')
legend('boxoff')


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

% Plot the correlations within a muscle 
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
filename2 = ['results_',datestr(now,'yyyymmddTHHMMSS'),'.mat'];
filename3 = ['sschain_',datestr(now,'yyyymmddTHHMMSS'),'.mat'];

save(filename,'chain');
save(filename2,'results'); 
save(filename3,'sschain')

user = memory;
memAfter = user.MemUsedMATLAB;

% display some stuff about the MCMC run
disp(['========================== ' ])
disp(['elapsed time = ' num2str(runtime/60) ' min'])
disp(['change in Matlab memory use = ' num2str((memAfter - memBefore)/1e6) ' MB'])
disp('   ')

%%
% Then re-run starting from the results of the previous run,
% this will take couple of minutes.
% options.nsimu = 10;
% [results, chain, s2chain] = mcmcrun(model,data,params,options, results);

function ss = Arm16_SS(theta,data)
 
%     time   = data.xdata;
%     ydata  = data.ydata(:,2:end);
    ydata_pos = data.ydata(:,2); % this is for only caring about position
    ydata_vel = data.ydata(:,3);
    % we know initial state is pos = 1, vel = 0, but could add this as
    % variability later... 
       
    y = Arm16_fun(theta(1:60));
    pre_temp = theta(61); 
    
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


      offset = 4; % try between 2 and 4
      amplitude = 0.4; % try between 0 and 1
    
%     ss = likelihood + prior; 
    denominator = temperature(pre_temp,offset,amplitude);
    ss = (sum(likelihood))/denominator;  

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

function prior_sum = prior_act(theta,mu,sig)
    tfinal = 0.5;
    h = 0.001;
    time = (0:h:tfinal)';
    Muscle_1 = theta(1:10);
    Muscle_2 = theta(11:20);
    Muscle_3 = theta(21:30);
    Muscle_4 = theta(31:40);
    Muscle_5 = theta(41:50);
    Muscle_6 = theta(51:60);
    
    pre_temp = theta(61); 
    
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

    prior1 = ((SumIntegExc - 0.0)/0.08)^2;
     
    
    prior2 = (pre_temp)^2;
    
    
    offset = 4; % try between 2 and 4
    amplitude = 0.4; % try between 0 and 1
    denominator = temperature(pre_temp,offset,amplitude);
    
    prior_sum = (prior1/denominator) + prior2; 
    
    % post calc when temp is between less than 1.001. 
end 

function data = temperature(pre_temp,offset,amplitude) 
    
    data = 1 + inverse_logit(pre_temp-offset) * amplitude; 
    
end

function output = inverse_logit(input)
    if input>0 
        output = 1./(1+exp(-input));
    else
        output = exp(input)./(1+exp(input));
    end
end

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