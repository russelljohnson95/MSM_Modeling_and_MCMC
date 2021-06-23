function  rank_plot_SMD(chain,options,n_pools)
% rank plots 

% called by Spring_mass_Damper_mcmc_parallel.m 

% Author: Russell T. Johnson, rtjohnso@usc.edu
% Last Edited: 6-18-21

% ---------------------------------------------------

chain_one = chain(:,:,1);
chain_two = chain(:,:,2);
chain_thr = chain(:,:,3); 
chain_for = chain(:,:,4);
chain_fiv = chain(:,:,5); 

for j = 1:7 
   m_result = vertcat(chain_one(:,1),chain_two(:,1), chain_thr(:,1),chain_for(:,1), chain_fiv(:,1));
   c_result = vertcat(chain_one(:,2),chain_two(:,2), chain_thr(:,2),chain_for(:,2), chain_fiv(:,2));
   k1_result = vertcat(chain_one(:,3),chain_two(:,3), chain_thr(:,3),chain_for(:,3), chain_fiv(:,3));
   k2_result = vertcat(chain_one(:,4),chain_two(:,4), chain_thr(:,4),chain_for(:,4), chain_fiv(:,4));
   R_result = vertcat(chain_one(:,5),chain_two(:,5), chain_thr(:,5),chain_for(:,5), chain_fiv(:,5));
   x0_result = vertcat(chain_one(:,6),chain_two(:,6), chain_thr(:,6),chain_for(:,6), chain_fiv(:,6));
   xdot0_result = vertcat(chain_one(:,7),chain_two(:,7), chain_thr(:,7),chain_for(:,7), chain_fiv(:,7));
end

% nsimu = options.nsimu;
nsimu = size(chain,1);
Run = ones(nsimu*n_pools,1);
for i = 2:n_pools
    Run(((i*1-1)*nsimu+1):nsimu*i) = i;
end

m_result = ([Run,m_result]);
c_result = ([Run,c_result]);
k1_result = ([Run,k1_result]);
k2_result = ([Run,k2_result]);
R_result = ([Run,R_result]);
x0_result = ([Run,x0_result]);
xdot0_result = ([Run,xdot0_result]);

m_result_sort = sortrows(m_result,2);
c_result_sort = sortrows(c_result,2);
k1_result_sort = sortrows(k1_result,2);
k2_result_sort = sortrows(k2_result,2);
R_result_sort = sortrows(R_result,2);
x0_result_sort = sortrows(x0_result,2);
xdot0_result_sort = sortrows(xdot0_result,2);

Rank = (1:nsimu*n_pools)';
Rank_Norm = Rank./(nsimu*n_pools);

m_results_rank = ([m_result_sort,Rank_Norm]);
c_results_rank = ([c_result_sort,Rank_Norm]);
k1_results_rank = ([k1_result_sort,Rank_Norm]);
k2_results_rank = ([k2_result_sort,Rank_Norm]);
R_results_rank = ([R_result_sort,Rank_Norm]);
x0_results_rank = ([x0_result_sort,Rank_Norm]);
xdot0_results_rank = ([xdot0_result_sort,Rank_Norm]);

M_rank = sortrows(m_results_rank,1);
C_rank = sortrows(c_results_rank,1);
K1_rank = sortrows(k1_results_rank,1);
K2_rank = sortrows(k2_results_rank,1);
R_rank = sortrows(R_results_rank,1);
X0_rank = sortrows(x0_results_rank,1);
XDOT0_rank = sortrows(xdot0_results_rank,1);

figure()

subplot(2,4,1)
histogram(M_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(M_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(M_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(M_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(M_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Mass')

subplot(2,4,2)
histogram(C_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(C_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(C_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(C_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(C_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Damping')

subplot(2,4,3)
histogram(K1_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(K1_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(K1_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(K1_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(K1_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Stiffness 1')

subplot(2,4,4)
histogram(K2_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(K2_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(K2_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(K2_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(K2_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Stiffness 2')

subplot(2,4,5)
histogram(R_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(R_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(R_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(R_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(R_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Spring Transition')

subplot(2,4,6)
histogram(X0_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(X0_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(X0_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(X0_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(X0_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Init Pos')

subplot(2,4,7)
histogram(XDOT0_rank(1:nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
hold on 
histogram(XDOT0_rank(nsimu+1:2*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
hold on 
histogram(XDOT0_rank(2*nsimu+1:3*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
hold on 
histogram(XDOT0_rank(3*nsimu+1:4*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','g','FaceAlpha',0.5)
hold on 
histogram(XDOT0_rank(4*nsimu+1:5*nsimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','c','FaceAlpha',0.5)
title('Init Vel')

end

