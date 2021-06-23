% Make results plots for the Paper. Calculated from Results from an MCMC run of
% the mass-spring-damper system

% Author: Russell T. Johnson, rtjohnso@usc.edu
% Last Edited: 6-18-21

% ---------------------------------------------------

clear
clc
close all

% load results from latest chain
load chain_20210526T111546 % 30k iterations, 45 mins
load results_20210526T111546

% specify random number
rng(99)

% chain is the variable that is saved
titles = ["m", "c", "k1", "k2", "T", "x0", "xdot0"];
n_pools = size(chain,3);
choose_result = 2;  % specify which of the runs we want to pull data from 
n_iter = size(chain,1); 
burn_in = n_iter*0.50; % specify the burn in 
result_1 = results(:,:,1);
result_2 = results(:,:,2);
result_3 = results(:,:,3);
result_4 = results(:,:,4);
result_5 = results(:,:,5);

% read in prior
prior_center(1,:) = result_1.prior(:,1);
prior_width(1,:)  = result_1.prior(:,2); 

prior_center(2,:) = result_2.prior(:,1);
prior_width(2,:)  = result_2.prior(:,2); 

prior_center(3,:) = result_3.prior(:,1);
prior_width(3,:)  = result_3.prior(:,2); 

prior_center(4,:) = result_4.prior(:,1);
prior_width(4,:)  = result_4.prior(:,2); 

prior_center(5,:) = result_5.prior(:,1);
prior_width(5,:)  = result_5.prior(:,2); 

limits = result_1.limits;

% combine the five (or n_pools #) chains for results

non_burnin = n_iter - burn_in; 
chain_full = zeros(non_burnin*n_pools,size(chain,2)); 
for i = 1:n_pools 
    chain_full(1+(non_burnin*(i-1)):non_burnin*i,:) = chain(burn_in+1:end,:,i); 
    
end

% plot correlations 
figure(); clf
set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 14])
mcmcplot(chain(burn_in+1:end,:,choose_result),[],result_1,'pairs');

options.titles = titles; 
options.nsimu = size(chain,1); 

nsimu = size(chain,1) - burn_in;

real = [1, 0.20, 3, 10, 0.03, 0.1, 0];

% --------- Here's the main figure ------------------- % 

figure(99)

set(gcf,'units','centimeters','Position',[7.5935 4.2863 20 14])

ylimmin2 = [0, 0,  0,  0, 0, 0, -1];
ylimmax2 = [2, 1, 10, 20, .1, .2,  1];

% second row is the chain series from a single chain
for i = 1:7
    hAx(i) = subplot(5,8,8+i);
    for j = 1:n_pools
        h1 = plot(chain(:,i,1),'color','k','LineWidth',0.1);
        hold on 
        h2 = plot(chain(:,i,2),'color','b','LineWidth',0.1);
        hold on 
        h3 = plot(chain(:,i,3),'color','r','LineWidth',0.1);
        hold on 
        h4 = plot(chain(:,i,4),'color','g','LineWidth',0.1);
        hold on 
        h5 = plot(chain(:,i,5),'color','c','LineWidth',0.1);
        hold on 
    end
    h6 = xline(burn_in,'LineStyle','--','color',[ 0.4660    0.6740    0.1880],'LineWidth',1.5);
    set(hAx(i),'xtick',[])
    set(hAx(i),'xticklabel',[])
    yMaxValue = max(chain(:,i,choose_result))*1.2; % give the density some space to breathe
    ylim([ ylimmin2(i) ylimmax2(i)]); 
    if i == 1
        ylabel('value')
        text(-0.9,1,'B','fontsize',10,'fontweight','bold','units','normalized')
    end
    if i == 5
        yticks([0 0.1])
    end
    if i == 7
        lgd = legend([h1, h2, h3, h4, h5, h6],'Chain 1','Chain 2','Chain 3','Chain 4','Chain 5','burn-in','orientation','vertical','position',[0.82 0.6 0.1124 0.1597]); %width = 0.1597
        legend('boxoff')
        lgd.FontSize = 8;
    end
    set(hAx(i),'fontsize',9)
    box off

end


xlimmin3 = [0, 0,  0,  0, 0, -.10, -1];
xlimmax3 = [3, 1, 8, 20, 0.1, 0.4,  1];
% Third row is for the posterior denisty values --
for i = 1:7
    subplot(5,8,16+i)

    [y(:,i),x(:,i)]=density(chain_full(:,i),[]);
    plot(x(:,i),y(:,i),'-k','LineWidth',2)
    set(gca,'fontsize',9)
    xline(real(i),'r--','LineWidth',1.5);
    yMaxValue = max(y(:,i))*1.2; % give the density some space to breathe
    ylim([0 ceil(yMaxValue)])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    xlim([xlimmin3(i) xlimmax3(i)]);
    if i == 1
        ylabel('density')
        text(-0.9,1,'C','fontsize',10,'fontweight','bold','units','normalized')
    end
    if i == 7
        lgd= legend('Post. Density','Real Value','Orientation','vertical','position',[0.82 0.55 0.1124 0.07]);
        legend('boxoff')
        lgd.FontSize= 8;
    end
    box off
end 

% First row is the prior density distribution, with the real value and
% initial value for each parameter

% find the max values in the densities above - in other words - find the
% peak probability 
for i = 1:7
    y_max = max(y(:,i));
    B(:,i) = y_max == y(:,i);
    index(i) = find([B(:,i)] == 1);
    posterior_center(i) = x(index(i),i); 
end

% Need to find the standard deviation - which i think this is the best way
% to do that... It's kind of wonky - but the idea is to integrate the AUC
% of the posterior density, until we get to 34.1%. This sort of assumes
% that it's a normal distribution, which doesn't hold up exactly. 

integrate = zeros(100,7);
integrate_sum = zeros(100,7);

for i = 1:7 
    for j = index(i):100
        integrate(j-index(i)+1,i) = (x(j,i)-x(j-1,i)) * y(j,i);
    end
    integrate_sum(1,i) = integrate(1,i);
    for k = 2:100
        integrate_sum(k,i) = integrate_sum(k-1,i)+integrate(k,i);
    end
    C(:,i) = integrate_sum(:,i)<0.341;
    width(i) = index(i) + sum(C(:,i));
    posterior_std(i) = x(width(i),i)-posterior_center(i);
end  

xlimit_get = [xlimmin3; xlimmax3];
init_param = chain(1,:,choose_result);


for i = 1:7
    subplot(5,8,i)
%     ii  = inds(i);
%     mus = prior(:,i);
    mu  = prior_center(1,i);
    sig = prior_width(1,i);
    mi  = limits(i,1);
    ma  = limits(i,2);

    xlimit=xlimit_get(:,i);
    % xp = linspace(max(mi,mu-3*sig),min([ma,mu+3*sig,x(end)+xdr]));
    xp = linspace(max(mi,mu-3*sig),min([ma,mu+3*sig]));
    yp = norpf(xp,mu,sig^2);
    yn = nordf((mi-mu)/sig)+1-nordf((ma-mu)/sig); % area outside bounds
    plot([xp(1),xp],[0,yp./(1-yn)],'--k','LineWidth',1.5)
    xline(init_param(i),'color',[0.2422    0.1504    0.6603],'LineWidth',1.5);
    hold on 
    xline(real(i),'r--','LineWidth',1.5);
    ymaxval = max(yp./(1-yn))*1.2;
    ylim([0 ceil(ymaxval)]);
%     yticks([0 ceil(ymaxval)]);
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'fontsize',9)
    hold off
    xlim([xlimmin3(i) xlimmax3(i)]);
    title(titles{i})
    if i == 1
        ylabel('density')
        text(-0.9,1,'A','fontsize',10,'fontweight','bold','units','normalized')
    end
    if i == 4
        xticks([0 20])
    end
    if i == 7 
        lgd= legend('Prior Density','Init. Param.','Real Value','Orientation','vertical','position',[0.82 0.85 0.1124 0.10]);
        legend('boxoff')
        lgd.FontSize= 8;
    end
    box off
end

% Forth row is to show convergence across 5 chains. 
for i = 1:7
    %calculate the rank for each parameter 
    rank = rank_plot_MSD_param(chain(burn_in+1:end,i,:),n_pools);
      
    subplot(5,8,24+i)
    set(gca,'fontsize',16)
    histogram(rank(1:nsimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
    hold on 
    histogram(rank(nsimu+1:2*nsimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
    hold on 
    histogram(rank(2*nsimu+1:3*nsimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
    hold on 
    histogram(rank(3*nsimu+1:4*nsimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
    hold on 
    histogram(rank(4*nsimu+1:5*nsimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
    set(gca,'fontsize',9)
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    box off
    
    if i == 1
        text(-0.9,1,'D','fontsize',10,'fontweight','bold','units','normalized')
    end
    if i == 7
        lgd = legend('Run 1','Run 2','Run 3','Run 4','Run 5','orientation','vertical','position',[0.82 0.25 0.1124 0.1597]);
        legend('boxoff')
        lgd.FontSize= 8;
    end
end

% do 10 random draws from each chain, then plot the results to see how
% the results fit the data
time = 0:0.05:20; % time column, !! could make more flexible based on input data
% Load Data from results
load('position_newmag20.mat')
load('position_4trials_plusnoise.mat')
for i = 1:size(position_noise,2)
    velocity_noise(:,i) = deriv(position_noise(:,i),0.05);
end

for k = 1:20
   draw(k) = randi([burn_in+1 options.nsimu]);
   Draw_Results(k,:) = chain(draw(k),:,choose_result);
%        y0 = [2,0];
%      y0 = [Draw_Results(k,4),0];
   y0(k,:) = [Draw_Results(k,6),Draw_Results(k,7)];

   [t(:,k),oscillator(:,k*2-1:k*2)] = ode15s(@MSD_sys,time,y0(k,:),[],Draw_Results(k,:));
end

% last row -  plot the kinematics from the random draws 

subplot(5,8,35:37)
for j = 1:19
    plot(time,oscillator(:,j*2-1),'color',[.17 .17 .17],'LineWidth',2)
    hold on 
end
h2 = plot(time,oscillator(:,20*2-1),'color',[.17 .17 .17],'LineWidth',2);
hold on
h1 = plot(time,position,'g','LineWidth',1.5);
h3 = plot(time,position_noise(:,1),'r--','LineWidth',2);
hold on
for k = 2:4
     plot(time,position_noise(:,k),'r--','LineWidth',2);
end
lgd = legend([h1 h2 h3],'Ref. Position','Random Draws','Pos. with Noise','Orientation','vertical');
legend('boxoff')
lgd.FontSize= 8;
box off
set(gca,'fontsize',9)
text(-0.888,1,'E','fontsize',10,'fontweight','bold','units','normalized')
ylim([-0.2 0.2])
xlabel('time (s)')
ylabel('position (m)')


% mechanics functions to calculate the motion of the mass
function y=MSD_fun(time,theta,y0)


    [t,y] = ode15s(@MSD_sys,time,y0,[],theta);

end

function ydot = MSD_sys(t,y,theta)
    % ode system function for MCMC mass spring damper example
    
    m = theta(1);
    c = theta(2); 
    k1 = theta(3); % stiffness value for region 1
    k2 = theta(4); % stiffness value for region 2
    R  = theta(5);

    
    if (-R<y(1)) && (y(1)<R)
        A = [0 1; -k1/m -c/m];
    else 
        A = [0 1; -(((k1*R)/(abs(y(1)))-(k2*R)/(abs(y(1)))+k2))/m -c/m];

    end
    
    ydot = A*y ;

end



