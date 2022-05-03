function rank_plot_Arm16_10CRBVs_fxn(chain)

% load('chain_results_20200814T094133.mat');

muscles = {'Tri Long', 'Tri Lat','Tri Med','Biceps LH', 'Biceps SH','Brachior'};

nMuscles = 6; 
nNodes = 10; 

for i = 1:nMuscles

    rank_plot_muscle_10node(chain(:,((i*nNodes)-9):(i*nNodes),:),muscles{i})

end

% nSimu = size(chain,1);
% nParams = size(chain,2);
% nChains = size(chain,3);
% Rank = (1:nSimu*nChains)';
% Rank_Norm = Rank./(nSimu*nChains);
% 
% Run = ones(nSimu*nChains,1);
% for i = 2:nChains
%     Run(((i*1-1)*nSimu+1):nSimu*i) = i;
% end

% % position 
% pos = vertcat(chain(:,49,1),chain(:,49,2), chain(:,49,3));
% pos_run = ([Run,pos]);
% pos_run_sort = sortrows(pos_run,2);
% pos_rank = ([pos_run_sort,Rank_Norm]);
% pos_rank_out = sortrows(pos_rank,1);
% 
% figure()
% histogram(pos_rank_out(1:nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
% hold on 
% histogram(pos_rank_out(nSimu+1:2*nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
% hold on 
% histogram(pos_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
% title('position')
% 
% % velocity 
% vel = vertcat(chain(:,50,1),chain(:,50,2), chain(:,50,3));
% vel_run = ([Run,vel]);
% vel_run_sort = sortrows(vel_run,2);
% vel_rank = ([vel_run_sort,Rank_Norm]);
% vel_rank_out = sortrows(vel_rank,1);
% 
% figure()
% histogram(vel_rank_out(1:nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','k','FaceAlpha',0.5)
% hold on 
% histogram(vel_rank_out(nSimu+1:2*nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','b','FaceAlpha',0.5)
% hold on 
% histogram(vel_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',20,'BinLimits',[0 1],'FaceColor','r','FaceAlpha',0.5)
% title('velocity')




end
