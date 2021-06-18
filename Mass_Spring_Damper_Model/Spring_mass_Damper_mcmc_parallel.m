clear
clc
close all

user = memory;
memBefore = user.MemUsedMATLAB;

rng(99)

% run a mcmc algorithm to predict/recover the mass, stiffness, and damping values
% for a mass-spring-damper system based on some data

% written by Russell T Johnson, University of Southern California
% rtjohnso@usc.edu
% Last edited: 5/18/21

% This depends on following MCMC package for MATLAB: 
% https://mjlaine.github.io/mcmcstat/

% True Model Values
m_real  = 1; % mass
c_real  = 0.20; % damping 
k1_real = 3; % stiffness value for region 1
k2_real = 10; % stiffness value for region 2
T_real  = 0.03; % threshold to switch from region 1 to region 2
x0_real = 0.1; % inital position
xdot0_real = 0; % initial velocity 

% Load Data from results
load('position_4trials_plusnoise.mat')

% find derivative 
velocity_noise = zeros(size(position_noise,1),size(position_noise,2)); 
for k = 1:size(position_noise,2)
    velocity_noise(:,k) = deriv(position_noise(:,k),0.05);
end

% define time column for simulation
time = 0:0.05:20; 

%% set up data

data.xdata = time';
data.ydata = [time',position_noise,velocity_noise];

%% Create Model

model.ssfun = @SMD_SS;

%% this is the part where we do the parallelization
% must have parallel tool box installed in matlab - check first!

check = contains(struct2array(ver), 'Parallel Computing Toolbox');

if check ~=1
    error('check to see parallel toolbox is installed in matlab')
end
  
% set number of parallel pools 
n_pools = 5; 

% due to simulation: create some variances based on assumptions for
% measurement accuracy - this gets used to set priors differently for each
% of the 5 parallel chains 
m_var     = 0.30; 
c_var     = 0.08;
k1_var    = 0.50; 
k2_var    = 1.00; 
T_var     = 0.02; 
x0_var    = 0.05; 
xdot0_var = 0.20; 

% this adds 10% noise to the prior center
prior_center = (randn(7,n_pools));
prior_noise  = 0.1; 

for j = 1:n_pools
    m_prior(j) = m_real + (prior_center(1,j) * (m_var*prior_noise));
    c_prior(j) = c_real + (prior_center(2,j) * (c_var*prior_noise));
    k1_prior(j) = k1_real + (prior_center(3,j) * (k1_var*prior_noise));
    k2_prior(j) = k2_real + (prior_center(4,j) * (k2_var*prior_noise));
    T_prior(j) = T_real + (prior_center(5,j) * (T_var*prior_noise));
    x0_prior(j) = x0_real + (prior_center(6,j) * (x0_var*prior_noise));
    xdot0_prior(j) = xdot0_real  + (prior_center(7,j) * (xdot0_var*prior_noise)); % can't multiply by 0
end

% define initial parameters based on prior for each chain
for j = 1:n_pools
    m(j)     = m_prior(j) + randn(1,1) * m_var;
    c(j)     = c_prior(j) + randn(1,1) * c_var;
    k1(j)    = k1_prior(j) + randn(1,1) * k1_var;
    k2(j)    = k2_prior(j) + randn(1,1) * k2_var;
    T(j)     = T_prior(j) + randn(1,1) * T_var;
    x0(j)    = x0_prior(j) + randn(1,1) * x0_var;
    xdot0(j) = xdot0_prior(j) + randn(1,1) * xdot0_var;
end

% if any of model parameter numbers are less than zero, the dynamics of the system
% breaks and everything dies, so check to see if the initial proposals come
% out to be < 0 and then reset them to arbitrary positive values. 

for j = 1:n_pools
    if m(j) <0
        m(j)  = m_real  + randn(1,1)* 0.05;
    end
    if c(j) < 0
        c(j)  = c_real  + randn(1,1)* 0.05;
    end
    if k1(j) < 0
        k1(j) = k1_real + randn(1,1)* 0.1;
    end
    if k2(j) < 0 
        k2(j) = k2_real + randn(1,1)* 0.1; 
    end
    if T(j) < 0
        T(j)  = T_real  + randn(1,1)* 0.1;
    end
end

titles = ["m", "c", "k1", "k2", "T", "x0", "xdot0"];

for i = 1:n_pools
    params(:,i) = { 
        {'m',    m(i),       0, 50,  m_prior(i),     (m_var*2)} 
        {'c',    c(i),       0, 10,  c_prior(i),     (c_var*2)}
        {'k1',   k1(i),      0, 50,  k1_prior(i),    (k1_var*2)}
        {'k2',   k2(i),      0, 50,  k2_prior(i),    (k2_var*2)}
        {'T',    T(i),       0,  2,  T_prior(i),     (T_var*2)}
        {'x0',   x0(i)    -2,  2,  x0_prior(i),    (x0_var*2)}
        {'xdot0',xdot0(i),-20, 20, xdot0_prior(i), (xdot0_var*2)}
        };    
end

% set up the options for mcmcrun 
options.nsimu = 30000;
options.stats = 1; 
options.stats2 = 1; 
options.waitbar = 0;

% define burn-in as a percentage of the number of simulations... 
burn_in = options.nsimu *0.50;

% start the clock
tic

% Open the parallel pools 
poolobj = parpool(n_pools);

% ----------------- DO THE MCMC!!!! -------------------------------
parfor k = 1:n_pools
    [results(:,:,k), chain(:,:,k), s2chain(:,:,k)]= mcmcrun(model,data,params(:,k),options);
end

% Stop the clock
runtime = toc;

% stop the parallel pools
delete(poolobj)

% assemble the results
results_one = results(:,:,1);
results_two = results(:,:,2);
results_thr = results(:,:,3); 
results_for = results(:,:,4);
results_fiv = results(:,:,5); 

% plot the chains
figure(9);
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 12])
for i =1:size(chain,2)
    subplot(2,4,i)
    plot(chain(:,i,1),'k')
    hold on 
    plot(chain(:,i,2),'b')
    hold on 
    plot(chain(:,i,3),'r')
    hold on 
    plot(chain(:,i,4),'g')
    hold on 
    plot(chain(:,i,5),'c')
    title(titles{i})
    xlabel('iteration')
    ylabel('value')
    box off
end

% plot more results from MCMC toolbox - this creates a bunch of plots,
% uncomment if you want to look at each of the chains
% for i = 1:size(chain,3)
%     figure(); clf
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])
%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'pairs');
% 
%     figure(); clf
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])
%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel',2);
%     hold on 
% 
%     figure(); clf
%     set(gcf,'units','centimeters','Position',[7.5935 4.2863 10 13])
%     mcmcplot(chain(burn_in:end,:,i),[],results(:,:,i),'denspanel2',2);
% end

% do 10 random draws from each chain, then plot the results to see how
% the results fit the data
for i = 1:n_pools
    for k = 1:10
       draw(k) = randi([burn_in options.nsimu]);
       Draw_Results(k,:,i) = chain(draw(k),:,i);
       y0(k,:,i) = [Draw_Results(k,6,i),Draw_Results(k,7,i)];
       [t(:,k),oscillator(:,k*2-1:k*2,i)] = ode15s(@SMD_sys,time,y0(k,:,i),[],Draw_Results(k,:,i));
    end
end

% plot the random draws 
figure()
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 13])
for i = 1:n_pools 
    subplot(3,2,i)
    for j = 1:10
        plot(time,oscillator(:,j*2-1,i),'color',[.17 .17 .17],'LineWidth',2)
        hold on 
    end
    h1 = plot(time,position_noise(:,1),'r--','LineWidth',2);
    hold on
    for k = 2:4
         plot(time,position_noise(:,k),'r--','LineWidth',2);
    end
    if i == 1
        legend([h1],'Ref. Positions')
        legend('boxoff')
    end
    set(gca,'fontsize',10)
    xlabel('time (s)')
    ylabel('position (m)')
    box off
end 

% Plot the rank without the first burn in 
rank_plot_SMD(chain(burn_in:end,:,:),options,n_pools);

% save the chains so that we can do other analysss later if we want
filename = ['chain_',datestr(now,'yyyymmddTHHMMSS'),'.mat'];
save(filename,'chain');
filename = ['results_',datestr(now,'yyyymmddTHHMMSS'),'.mat'];
save(filename,'results');

user = memory;
memAfter = user.MemUsedMATLAB;

% return

% display some stuff about the MCMC run
disp(['========================== ' ])
disp(['elapsed time = ' num2str(runtime/60) ' min'])
disp(['change in Matlab memory use = ' num2str((memAfter - memBefore)/1e6) ' MB'])
disp('   ')

%% End of main function, begin other functions

function ss = SMD_SS(theta,data)
    % sum of squares function for the posterior probability 

    time   = data.xdata;
    ydata  = data.ydata(:,2:end); % just takes the position data

    y0 = [theta(6),theta(7)];
    
    ymodel = SMD_fun(time,theta,y0);
    n_trials = 4; 
    
    for i = 1:n_trials
        sumsquare(i) = sum(((ymodel(:,1) - ydata(:,i))/0.10).^2);
    end
    ss = sum(sumsquare);
end

function y=SMD_fun(time,theta,y0)
    % integration of the mass spring damper system 
    [t,y] = ode15s(@SMD_sys,time,y0,[],theta);

end

function ydot = SMD_sys(t,y,theta)
    % ode system function for MCMC mass spring damper example
    
    m = theta(1);
    c = theta(2); 
    k1 = theta(3); % stiffness value for region 1
    k2 = theta(4); % stiffness value for region 2
    T  = theta(5);
    
    if (-T<y(1)) && (y(1)<T)
        A = [0 1; -k1/m -c/m];
    else 
        A = [0 1; -(((k1*T)/(abs(y(1)))-(k2*T)/(abs(y(1)))+k2))/m -c/m];
    end
    
    ydot = A*y ;

end
