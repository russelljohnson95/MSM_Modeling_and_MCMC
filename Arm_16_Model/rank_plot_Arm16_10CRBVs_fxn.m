function rank_plot_Arm16_10CRBVs_fxn(chain)

% written by Russell T Johnson, University of Southern California
% rtjohnso@usc.edu
% Last edited: 6/18/2021    

% more or less a wrapper function to send out to plot the ranks for each of
% the muscles in the model 

muscles = {'Tri Long', 'Tri Lat','Tri Med','Biceps LH', 'Biceps SH','Brachior'};

nMuscles = 6; 
nNodes = 10; 

for i = 1:nMuscles

    rank_plot_muscle_10node(chain(:,((i*nNodes)-9):(i*nNodes),:),muscles{i})

end

end
