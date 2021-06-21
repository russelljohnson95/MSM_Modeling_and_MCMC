function  rank = rank_plot_SMD_param(chain,n_pools)
% rank plots but for one parameter at a time

% called by Figure_Manuscript_fromChain.m

% Author: Russell T. Johnson, rtjohnso@usc.edu
% Last Edited: 6-18-21

% ---------------------------------------------------

chain_one = chain(:,:,1);
chain_two = chain(:,:,2);
chain_thr = chain(:,:,3); 
chain_for = chain(:,:,4);
chain_fiv = chain(:,:,5); 

result = vertcat(chain_one(:,1),chain_two(:,1), chain_thr(:,1),chain_for(:,1), chain_fiv(:,1));

% nsimu = options.nsimu;
nsimu = size(chain,1);
Run = ones(nsimu*n_pools,1);
for i = 2:n_pools
    Run(((i*1-1)*nsimu+1):nsimu*i) = i;
end

result = ([Run,result]);

result_sort = sortrows(result,2);

Rank = (1:nsimu*n_pools)';
Rank_Norm = Rank./(nsimu*n_pools);

results_rank = ([result_sort,Rank_Norm]);

rank = sortrows(results_rank,1);


end

