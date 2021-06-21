function [y] = Arm16_SimManager_controls_CRBF_6musc(theta)
% written by Russell T Johnson, University of Southern California
% rtjohnso@usc.edu
% Last edited: 6/18/2021

% This code takes in the theta (or amplitudes for the CRBFs), computes the
% corresponding muscle excitations, then uses OpenSim integration to
% compute the kinematics. 

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model
osimModel = Model('arm16_millard_rigidtendon.osim');

% muscle controller - this likely changes with OpenSim 4.x
muscleController = PrescribedController();
muscleController.setName('PiecewiseLinearFunction'); 
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

% amplitudes of each of the muscles 
Muscle_1 = theta(1:10);
Muscle_2 = theta(11:20);
Muscle_3 = theta(21:30);
Muscle_4 = theta(31:40);
Muscle_5 = theta(41:50);
Muscle_6 = theta(51:60);

Nrows = length(time);

controls = zeros(Nrows,Ncontrols);

% convert CRBF amplitudes to muscle excitations (controls)
controls(:,1) = CRBF_excit(time,Muscle_1);
controls(:,2) = CRBF_excit(time,Muscle_2);
controls(:,3) = CRBF_excit(time,Muscle_3);
controls(:,4) = CRBF_excit(time,Muscle_4);
controls(:,5) = CRBF_excit(time,Muscle_5);
controls(:,6) = CRBF_excit(time,Muscle_6);

% initial values of the system, in theory, this could be an input file or
% it could be a variable that gets used with the MCMC
initStateValues = zeros(8,1);

initStateValues(1,1) = 0;   % initial speed
initStateValues(2,1) = 0; % initial velocity
initStateValues(3,1) = 0.05; % initial activation tri long
initStateValues(4,1) = 0.05; % init act tri lat
initStateValues(5,1) = 0.05; % init act tri med
initStateValues(6,1) = 0.05; % init act biceps long head
initStateValues(7,1) = 0.05; % init act biceps short head 
initStateValues(8,1) = 0.05; % init act brachialis

numVar = osimState.getNY();
for i = 0:1:numVar-1
   osimState.updY().set(i, initStateValues(i+1,1));
end

osimModel.setPropertiesFromState(osimState);

% Get a reference to the controller that is attached to the model
muscleController = PrescribedController.safeDownCast(osimModel.getControllerSet.get(0));

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
osimModel.equilibrateMuscles(osimState);

% Record current state values before the simulation runs
numVar     = osimState.getNY();
initStates = zeros(numVar,1);
for i = 0:1:numVar-1
   initStates(i+1,1) = osimState.getY().get(i);
end

% Create a Manager to run the simulation
simulationManager = Manager(osimModel);
simulationManager.setInitialTime(0.0);
simulationManager.setFinalTime(tfinal);
simulationManager.setWriteToStorage(true);
%simulationManager.setPerformAnalyses(true);
simulationManager.setIntegratorAccuracy(1e-04)
simulationManager.integrate(osimState);

% Pull out the time and muscle activations
TimeArray = ArrayDouble();
StatesArray = ArrayDouble();

simulationManager.getStateStorage().getTimeColumn(TimeArray,-1);
ArrayLength = TimeArray.getSize(); 
Time_SM = zeros(ArrayLength,1);
for i = 1:ArrayLength
   Time_SM(i) = TimeArray.getitem(i-1);
end

Kinematics = zeros(ArrayLength,2);
for i = 1:2
   simulationManager.getStateStorage().getDataColumn(i-1,StatesArray);
   for j = 1:ArrayLength
      Kinematics(j,i) = StatesArray.getitem(j-1);
   end
end

% return the time column, the positions, and velocities of the elbow
y = [Time_SM,Kinematics]; 

end

