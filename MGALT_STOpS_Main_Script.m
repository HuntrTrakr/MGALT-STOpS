%% STOpS - Low-Thrust Trajectory Optimization Main Script 

% House cleaning
close all; clear; clc   % Clears up the workspace and command window

% Add other folders
addpath     Algorithms\                         ...
            Algorithms\Algorithm_Parameters\    ...
            Algorithms\Differential_Evolution\	...
            Algorithms\Genetic_Algorithm\   	...
            Algorithms\Monotonic_Basin\       	...
            Algorithms\Particle_Swarm\        	...
            Data\                               ...
            Data\JPL_Horizons\                  ...
            Error_Messages\                     ...
            Island_Model\                       ...
            Options\                            ...
          	Other_Functions\                    ...
            Solvers\                            ...
            Solvers\Cost_Functions\             ...
            Solvers\Direct\                     ...
            Solvers\Indirect\                   ...
            Solvers\Transfer_Conditions\        ...
            Thrust_Profiles\
        
% Name of folder to save results into
% *This is on your current path
OPT.save_folder = 'Optimization_Results';
        
% Do you want to use parallel processing?
% WARNING: This uses a lot of RAM. 
% Not recommended unless the system has 16GB+ RAM installed
OPT.parallel = true;



%% Departure/Transfer/Arrival Planet and Earliest/Latest Departure

% Section Description:
%{
-Each body is a string cooresponding to one of the 9 planets.
-Options are: 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 
  'Uranus', and 'Neptune'.
-The position and velocity of the departure body at the start time will 
  be used to create the inital conditions for the trajectory. The position 
  and velocity of the arrival body at the end time will be used to create 
  the end conditions for the trajectory.
-If the trajectory's goal is not to arrive at the arrival body, but just 
  to insert into its orbit, the arrival body's location will be ignored 
  and the end conditions become the position and velocity at some point 
  on its orbit.
-All entries are spelling and case sensitive


-Here the user is required to enter the earliest and latest departure date
  from the departure body. The earliest date is refered to as window1 and
  the latest departure date is refered to as window2.  The format for the
  date is [YEAR,MONTH,DAY] with each input (year, month, and day) being an
  integer.

-Time of flight parameters for the whole solver in days.
  It is important to specify the desired transfer time between each
  transfer segment, with the allowable tolerance.
    -Ex: On a mission from Earth->Mars->Jupiter, there are 2 flight
    segments which need a time specified for them
        tof_total = [300;900]
        tof_margin = [50,50;250,250]
    This tells the solver that a trajectory time of 300+-50 days is desired
    for the transfer from Earth->Mars, and 900+-250 days is desired from
    Mars->Jupiter
-Add margin to the desired days for random perturbations
%}


% ***User selection***
BOD.bodies = {'Earth', 'Mars', 'Jupiter'};
OPT.tof_total = [   700; ...
                    1500];       % tof w/ gravity assists [days]
OPT.tof_margin = [  200,200; ...
                    500,500];	% margin on tof -/+     [days]

% % Launch Windows
BOD.window1 = [1990,1,1];       % earliest launch day	[year,month,day]
BOD.window2 = [1995,1,1];       % latest launch day  	[year,month,day]



%% Cost Function Selection

% Section Description:
%{
-The chosen cost funtion dictates how the trajectory is described. A 
  variable string decribes each possible trajectory. In this work there 
  are two main structures for the variable strings.

-In what will be called the segmented method the trajectory is divided 
  into N segments each with a thrust value and a thrust pointing angle. 
  The variable string for the segmented method consists of a departure 
  time, the thrust and pointing angle for each segment, and an arrival 
  time. 

-In what will be called the costate method the trajectory has a consant 
  or equation thrust and the time dependent thrust angle is described 
  with three intial costate variables. The variable string for the 
  Conway method consists of a departure time, the three costate variables, 
  and an arrival time. 

-The two cost function handles available coorespond to the two variable 
  composition methods. They are: 'EP_cost_fun_segmented_2D' and 
  'EP_cost_fun_costate_2D'. The cost function is spelling and case 
  sensitive.
%}


% ***User selection***

% Solver selection
% OPT.solver = 'LT_DIR_FSM_2D';
% OPT.solver = 'LT_IN_FSM_2D';
% OPT.solver = 'MGALT_DIR_FBSM_2D';
OPT.solver = 'MGALT_IN_FBSM_2D';



%% Optimization Options

% Thrust Profile     
OPT.thrust = optionsLowThrust(OPT,'Yam-STOUR');
        
% Island Parameters
OPT.island = optionsIsland('3M-all');

% Cost and Weighting Parameters
[OPT.cost,OPT.weighting,OPT.ode] = optionsCost(OPT,'default');

% Algorithm Parameters
OPT.GA(1,1) = parametersGA('75_30');        % Genetic Algorithm Island #1
OPT.DE(1,1) = parametersDE('75_30');        % Differential Evolution Island #1
OPT.PSO(1,1) = parametersPSO('50_50');      % Particle Swarm Island #1
OPT.MBH(1,1).N1_Outer = 2000;
OPT.MBH(1,1).N1_Inner = 500;
OPT.MBH(1,1).N2_Outer = 250;
OPT.MBH(1,1).N2_Inner = 100;
OPT.MBH(1,1).maxclst = 5;
OPT.MBH(1,1).per_feas = [3.00,1.75,1.25,1.00];
OPT.MBH(1,1).per_rand = [0.1,0.07,0.04,0.01];
OPT.MBH(1,1).feas_check = @any;
OPT.MBH(1,1).feas_tol = 1e-5;



%% User Error Catch

% Check for errors
if errorMain(BOD,OPT)
    return
end


        
%% Pre Optimization

% Start parallel processing
if OPT.parallel
    if isempty(gcp())
        parpool('local')
    end
	OPT.ppool = gcp();
end

% Define the search window for the planetary vectors
BOD.windows = [BOD.window1; BOD.window2];
BOD.JD1 = getJulianDate(BOD.windows(1,1:3));	% earliest Julian date
BOD.JD2 = getJulianDate(BOD.windows(2,1:3)) + sum(OPT.thrust.tt_end) + ...
    sum(OPT.thrust.time(:,end));                % latest Julian Date
BOD.tspan = [BOD.JD1,BOD.JD2];

% Define Position, Velocity, and Julian Date for Bodies
[BOD.bodies_R,BOD.bodies_V,BOD.bodies_JD] = planetBodies(BOD);

% Create variable limits
VAR = MGALT_varLimits(BOD,OPT);

% Get constants
CONST = getConstants(BOD);



%% Optimization

% Section Description:
%{
-Here is where the optimization is actually run. All of the user inputs
  from the previous sections are run through here into the optimization
  algorithms.
%}

% Define selected as an empty struct
selected = struct;

% Timer
tic     % Start a timer

% Run through the solver
for num_mig = 1:OPT.island.Nmig+1
    
    % Algorithm counter in case an algorithm is used more than once
	count_GA = 1; 
    count_DE = 1; 
    count_PSO = 1; 
    count_MBH = 1;
    
    % Run through all the islands
	for num_isl = 1:OPT.island.Nisl
        
        switch char(OPT.island.isl_list(num_isl,:))
            
            case {'GA'}
                
                [eval_info(num_mig,num_isl), selected] = MGALT_GA(BOD,CONST,OPT,...
                    VAR,selected,num_mig,num_isl,count_GA); %#ok<*SAGROW>
                count_GA = count_GA + 1;
                
            case {'DE'}
                
                [eval_info(num_mig,num_isl), selected] = MGALT_DE(BOD,CONST,OPT,...
                    VAR,selected,num_mig,num_isl,count_DE); %#ok<*SAGROW>
                count_DE = count_DE + 1;
                
            case {'PSO'}
                
                [eval_info(num_mig,num_isl), selected] = MGALT_PSO(BOD,CONST,OPT,...
                    VAR,selected,num_mig,num_isl,count_PSO); %#ok<*SAGROW>
                count_PSO = count_PSO + 1;

            case {'MBH'}

                [eval_info(num_mig,num_isl), selected] = MGALT_MBH(BOD,CONST,OPT,...
                    VAR,selected,num_mig,num_isl,count_MBH); %#ok<*SAGROW>
                count_MBH = count_MBH + 1;
                
            otherwise
                
                errorPathDisplay()
                errorAlgorithm();
                return
                
        end
        
	end
    
end

run_time = toc/60;



%% Display Output

% Display all info to the Command Window and plot trajectory
displayResults(BOD,CONST,OPT,VAR,run_time,eval_info)



%% Save File

% Get current path with pwd and set new path
new_path = strcat(pwd,'\',OPT.save_folder);

% Get the solver name
switch OPT.solver
    
    case {'LT_IN_FSM_2D','LT_DIR_FSM_2D'}
        
        solver_name = OPT.solver(4:end-3);
        
    case {'MGALT_IN_FBSM_2D','MGALT_DIR_FBSM_2D'}
        
        solver_name = OPT.solver(7:end-3);
        
    otherwise
        
        errorPathDisplay()
        errorSolver()
        fprintf(2,'To save you data, manually save it from the workspace.\n');
        return
        
end

% Temp name
counter = 1;
temp_name = strcat(date(),'_',solver_name,'_',num2str(counter),'.mat');

% Check to see if name already exists
while isfile(fullfile(new_path,temp_name))
    counter = counter+1;
    temp_name = strcat(date(),'_',solver_name,'_',num2str(counter),'.mat');
end

% Save to the current path
save(temp_name)

% Move results to the name of "save_folder"
movefile(temp_name,new_path);

clear counter temp_name current_path new_path save_folder solver_name


