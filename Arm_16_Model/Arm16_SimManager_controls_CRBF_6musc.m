function [y] = Arm16_SimManager_controls_CRBF_6musc(theta)
% put some more description here if this all works 

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model

osimModel = Model('arm16_millard_rigidtendon.osim');
% osimModel = Model('arm16_millard_out.osim');

muscleController = PrescribedController();
muscleController.setName('PiecewiseLinearFunction'); % changed from 'PiecewiseLinear Controller' OS 3.3;
muscleController.setActuators(osimModel.updActuators());
osimModel.addController(muscleController);

% Initialize the model (this builds the system and initialize the state)
osimState = osimModel.initSystem();
Nstates    = osimModel.getNumStateVariables();
Ncontrols  = osimModel.getNumControls();
muscles    = osimModel.getMuscles(); 
nMuscles   = muscles.getSize();

tfinal = 0.5;
time = (0:0.005:tfinal)';
% time_star = -cos((k*pi)/(N-1));

Muscle_1 = theta(1:10);
Muscle_2 = theta(11:20);
Muscle_3 = theta(21:30);
Muscle_4 = theta(31:40);
Muscle_5 = theta(41:50);
Muscle_6 = theta(51:60);

Nrows = length(time);

controls = zeros(Nrows,Ncontrols);

controls(:,1) = CRBF_excit(time,Muscle_1);
controls(:,2) = CRBF_excit(time,Muscle_2);
controls(:,3) = CRBF_excit(time,Muscle_3);
controls(:,4) = CRBF_excit(time,Muscle_4);
controls(:,5) = CRBF_excit(time,Muscle_5);
controls(:,6) = CRBF_excit(time,Muscle_6);

% initStateValues = load('ArtificialData_Result_CRBF_6musc_states.txt')';
initStateValues = zeros(11,1);

% set initial state values sort of manually... 
% initStateValues(1,1) = 0;   % initial speed
% initStateValues(2,1) = 0; % initial velocity
% initStateValues(3,1) = 0.02; % init act tri long
% initStateValues(4,1) = 0.1366;
% initStateValues(5,1) = 0.02; % init act tri lat
% initStateValues(6,1) = 0.0705;
% initStateValues(7,1) = 0.02; % init act tri med
% initStateValues(8,1) = 0.0656;
% initStateValues(9,1) = 0.02; % init act biceps lh
% initStateValues(10,1) = 0.1487;
% initStateValues(11,1) = 0.02; % init act biceps sh
% initStateValues(12,1) = 0.1479;
% initStateValues(13,1) = 0.02; % init act brachioradialis 
% initStateValues(14,1) = 0.0934;

% set initial state values sort of manually... Then we can eventually work
% in the inital values to the MCMC algorithm

% With the rigid tendon muscle model model, we don't need to set fiber
% lengths

% initStateValues(1,1) = 0;   % initial speed
% initStateValues(2,1) = 0; % initial velocity
% initStateValues(3,1) = 0.05; % init act tri long
% initStateValues(4,1) = 0.1366;
% initStateValues(5,1) = 0.05; % init act tri lat
% initStateValues(6,1) = 0.0705;
% initStateValues(7,1) = 0.05; % init act tri med
% initStateValues(8,1) = 0.0656;
% initStateValues(9,1) = 0.05; % init act biceps lh
% initStateValues(10,1) = 0.1487;
% initStateValues(11,1) = 0.05; % init act biceps sh
% initStateValues(12,1) = 0.1479;
% initStateValues(13,1) = 0.05; % init act brachioradialis 
% initStateValues(14,1) = 0.0934;


initStateValues(1,1) = 0;   % initial speed
initStateValues(2,1) = 0; % initial velocity
initStateValues(3,1) = 0.02; % init act tri long
initStateValues(4,1) = 0.02;
initStateValues(5,1) = 0.02; 
initStateValues(6,1) = 0.02;
initStateValues(7,1) = 0.02; 
initStateValues(8,1) = 0.02;

numVar = osimState.getNY();
for i = 0:1:numVar-1
   osimState.updY().set(i, initStateValues(i+1,1));
end
% osimModel.equilibrateMuscles(osimState);
osimModel.setPropertiesFromState(osimState);

% Get a reference to the controller that is attached to the model
muscleController = PrescribedController.safeDownCast(osimModel.getControllerSet.get(0));

% Define functions for the control values for the muscles
% Tval = ArrayDouble(0,length(t));
% Eval = ArrayDouble(0,length(t));
numPoints = length(time);

Tval = ArrayDouble(0,numPoints);
Eval = ArrayDouble(0,numPoints);



for i = 1:nMuscles
   for j = 1:numPoints
      Tval.setitem(j-1,time(j,1));
      Eval.setitem(j-1,controls(j,i));
   end
   PLF = PiecewiseLinearFunction();
   for j = 1:numPoints
      PLF.addPoint(Tval.getitem(j-1),Eval.getitem(j-1));
   end
      muscleController.prescribeControlForActuator(i-1,PLF);
end



osimState = osimModel.initSystem();

% initActivations(1,1:nMuscles) = 0.02;
% 
% muscles = osimModel.getMuscles();
% for i = 0:nMuscles-1
%    muscles.get(i).setActivation(osimState, initActivations(i+1));
%    afl = ActivationFiberLengthMuscle.safeDownCast(muscles.get(i));
%    afl.setFiberLength(osimState, 0.1);
%    
% end
osimModel.equilibrateMuscles(osimState);

% Record current state values before the simulation runs
numVar     = osimState.getNY();
initStates = zeros(numVar,1);
for i = 0:1:numVar-1
   initStates(i+1,1) = osimState.getY().get(i);
end



% % Create a Manager to run the simulation
% simulationManager = Manager(osimModel);
% % simulationManager.setInitialTime(0.0);
% % simulationManager.setFinalTime(tfinal);
% simulationManager.setWriteToStorage(true);
% %simulationManager.setPerformAnalyses(true);
% % simulationManager.setIntegratorAccuracy(1e-04)
% simulationManager.integrate(tfinal);


% Create a Manager to run the simulation
simulationManager = Manager(osimModel);
simulationManager.setInitialTime(0.0);
simulationManager.setFinalTime(tfinal);
simulationManager.setWriteToStorage(true);
%simulationManager.setPerformAnalyses(true);
simulationManager.setIntegratorAccuracy(1e-04)
simulationManager.integrate(osimState);


% Uncomment to write out a states storage file
% simulationManager.getStateStorage().print('arm16_newnew_states.sto');


% Pull out the time and muscle activations
TimeArray = ArrayDouble();
StatesArray = ArrayDouble();

simulationManager.getStateStorage().getTimeColumn(TimeArray,-1);
ArrayLength = TimeArray.getSize(); 
Time_SM = zeros(ArrayLength,1);
for i = 1:ArrayLength
   Time_SM(i) = TimeArray.getitem(i-1);
end

% tic 

% States = zeros(ArrayLength,1);
% for i = 1:Nstates
%    simulationManager.getStateStorage().getDataColumn(i-1,StatesArray);
%    for j = 1:ArrayLength
%       States(j,i) = StatesArray.getitem(j-1);
%    end
% end


Kinematics = zeros(ArrayLength,2);
for i = 1:2
   simulationManager.getStateStorage().getDataColumn(i-1,StatesArray);
   for j = 1:ArrayLength
      Kinematics(j,i) = StatesArray.getitem(j-1);
   end
end

% position = zeros(ArrayLength,1);
% simulationManager.getStateStorage().getDataColumn(0,StatesArray);
% for j = 1:ArrayLength
%   position(j,1) = StatesArray.getitem(j-1);
% end

y = [Time_SM,Kinematics]; 

end

% toc
% 
% figure(2)
% subplot(2,1,1)
% plot(Time,rad2deg(States(:,1)));
% 
% subplot(2,1,2)
% plot(Time,rad2deg(States(:,2)));


