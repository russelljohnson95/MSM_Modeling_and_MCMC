function  rank_plot_muscle_10node(chain,title)
%UNTITLED3 Summary of this function goes here


nSimu = size(chain,1);
nChains = size(chain,3);

Run = ones(nSimu*nChains,1);
for i = 2:nChains
    Run(((i*1-1)*nSimu+1):nSimu*i) = i;
end

m1_p1 = vertcat(chain(:,1,1),chain(:,1,2), chain(:,1,3), chain(:,1,4), chain(:,1,5), chain(:,1,6), chain(:,1,7));
m1_p2 = vertcat(chain(:,2,1),chain(:,2,2), chain(:,2,3), chain(:,2,4), chain(:,2,5), chain(:,2,6), chain(:,2,7));
m1_p3 = vertcat(chain(:,3,1),chain(:,3,2), chain(:,3,3), chain(:,3,4), chain(:,3,5), chain(:,3,6), chain(:,3,7));
m1_p4 = vertcat(chain(:,4,1),chain(:,4,2), chain(:,4,3), chain(:,4,4), chain(:,4,5), chain(:,4,6), chain(:,4,7));
m1_p5 = vertcat(chain(:,5,1),chain(:,5,2), chain(:,5,3), chain(:,5,4), chain(:,5,5), chain(:,5,6), chain(:,5,7));
m1_p6 = vertcat(chain(:,6,1),chain(:,6,2), chain(:,6,3), chain(:,6,4), chain(:,6,5), chain(:,6,6), chain(:,6,7));
m1_p7 = vertcat(chain(:,7,1),chain(:,7,2), chain(:,7,3), chain(:,7,4), chain(:,7,5), chain(:,7,6), chain(:,7,7));
m1_p8 = vertcat(chain(:,8,1),chain(:,8,2), chain(:,8,3), chain(:,8,4), chain(:,8,5), chain(:,8,6), chain(:,8,7));
m1_p9 = vertcat(chain(:,9,1),chain(:,9,2), chain(:,9,3), chain(:,9,4), chain(:,9,5), chain(:,9,6), chain(:,9,7));
m1_p10 = vertcat(chain(:,10,1),chain(:,10,2), chain(:,10,3), chain(:,10,4), chain(:,10,5), chain(:,10,6), chain(:,10,7));

% Add column Run
m1_p1 = ([Run,m1_p1]);
m1_p2 = ([Run,m1_p2]);
m1_p3 = ([Run,m1_p3]);
m1_p4 = ([Run,m1_p4]);
m1_p5 = ([Run,m1_p5]);
m1_p6 = ([Run,m1_p6]);
m1_p7 = ([Run,m1_p7]);
m1_p8 = ([Run,m1_p8]);
m1_p9 = ([Run,m1_p9]);
m1_p10 = ([Run,m1_p10]);

% Sort from small to large
m1_p1_sort = sortrows(m1_p1,2);
m1_p2_sort = sortrows(m1_p2,2);
m1_p3_sort = sortrows(m1_p3,2);
m1_p4_sort = sortrows(m1_p4,2);
m1_p5_sort = sortrows(m1_p5,2);
m1_p6_sort = sortrows(m1_p6,2);
m1_p7_sort = sortrows(m1_p7,2);
m1_p8_sort = sortrows(m1_p8,2);
m1_p9_sort = sortrows(m1_p9,2);
m1_p10_sort = sortrows(m1_p10,2);

Rank = (1:nSimu*nChains)';
Rank_Norm = Rank./(nSimu*nChains);

% Add column labeling the rank
m1_p1_rank = ([m1_p1_sort,Rank_Norm]);
m1_p2_rank = ([m1_p2_sort,Rank_Norm]);
m1_p3_rank = ([m1_p3_sort,Rank_Norm]);
m1_p4_rank = ([m1_p4_sort,Rank_Norm]);
m1_p5_rank = ([m1_p5_sort,Rank_Norm]);
m1_p6_rank = ([m1_p6_sort,Rank_Norm]);
m1_p7_rank = ([m1_p7_sort,Rank_Norm]);
m1_p8_rank = ([m1_p8_sort,Rank_Norm]);
m1_p9_rank = ([m1_p9_sort,Rank_Norm]);
m1_p10_rank = ([m1_p10_sort,Rank_Norm]);

% sort
m1_p1_rank_out = sortrows(m1_p1_rank,1);
m1_p2_rank_out = sortrows(m1_p2_rank,1);
m1_p3_rank_out = sortrows(m1_p3_rank,1);
m1_p4_rank_out = sortrows(m1_p4_rank,1);
m1_p5_rank_out = sortrows(m1_p5_rank,1);
m1_p6_rank_out = sortrows(m1_p6_rank,1);
m1_p7_rank_out = sortrows(m1_p7_rank,1);
m1_p8_rank_out = sortrows(m1_p8_rank,1);
m1_p9_rank_out = sortrows(m1_p9_rank,1);
m1_p10_rank_out = sortrows(m1_p10_rank,1);

%plot
figure()

subplot(2,5,1)
histogram(m1_p1_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p1_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 1')

subplot(2,5,2)
histogram(m1_p2_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p2_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 2')

subplot(2,5,3)
histogram(m1_p3_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p3_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 3')

subplot(2,5,4)
histogram(m1_p4_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p4_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 4')

subplot(2,5,5)
histogram(m1_p5_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p5_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 5')

subplot(2,5,6)
histogram(m1_p6_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p6_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 6')

subplot(2,5,7)
histogram(m1_p7_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p7_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 7')

subplot(2,5,8)
histogram(m1_p8_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p8_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 8')

subplot(2,5,9)
histogram(m1_p9_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p9_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 9')

subplot(2,5,10)
histogram(m1_p10_rank_out(1:nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#0072BD','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(nSimu+1:2*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#D95319','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(2*nSimu+1:3*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#EDB120','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(3*nSimu+1:4*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#7E2F8E','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(4*nSimu+1:5*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#77AC30','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(5*nSimu+1:6*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#4DBEEE','FaceAlpha',0.5)
hold on 
histogram(m1_p10_rank_out(6*nSimu+1:7*nSimu,3),'NumBins',10,'BinLimits',[0 1],'FaceColor','#A2142F','FaceAlpha',0.5)
ylabel('Node 10')


sgtitle(title);

end

