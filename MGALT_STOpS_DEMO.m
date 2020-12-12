%% STOpS - MGALT Demo Read ME:

%{
This DEMO script is intended to give an introduction to using MGALT STOpS.

For this demo, the program will optimize a trajectory between Earth and
Jupiter with a Mars gravity assist using the indirect method.

Each section of the code will give a step by step instruction on what STOpS
is doing at that point, and will allow a first time user to understand the
code better.

User input is only required in the following sections:
"Initial"
"Departure/Transfer/Arrival Planet and Earliest/Latest Departure"
"Cost Function Selection"
"Optimization Options"
and the sections following that are all used by STOpS to perform error
checking, starting the optimization process, performing the optimization,
and saving the results.

At the end of the optimization process, the demo script will automatically
open another script used to load previous results.


*NOTE:
STOpS is not dependent on any toolboxes for MATLAB; however, it can
optionally use the "Parallel Computing in Optimization Toolbox" to speed up
the optimization process.
%}



%% Initial

% --- House cleaning:
% Closes open figures, clears the workspace from vars, and clears the
% command window
close all; clear; clc

% --- Add other folders:
% Adds the paths for the folders in "STOpS-MGALT" to the current file
% directory. This allows MATLAB to acces these folders, which are necesary
% for the optimization process
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

% --- Name of folder to save results into
% This is the name of the folder where results of the optimization process
% will be saved into. This folder will be on your current path, therefore,
% it will be contained as another folder in "STOpS-MGALT"
OPT.save_folder = 'Optimization_Results';

% --- Do you want to use parallel processing?
% WARNING: This uses a lot of RAM. 
% Not recommended unless the system has 16GB+ RAM installed
OPT.parallel = false;	% "true" ot "false"



%% Departure/Transfer/Arrival Planet and Earliest/Latest Departure

% --- Which planets do you want to optimize a trajectory between?
% The first planet in the cell, in this case "Earth", is the departure
% planet. The next planet, in this case "Mars", is the first gravity assist
% planet. The final planet, in this case "Jupiter", is the target planet. 
BOD.bodies = {'Earth', 'Mars', 'Jupiter'};

% --- Set the desired time of flight between each of the planets
% The first number, in this case "300" is the desired transfer time between
% Earth and Mars in days, and "1500" is the desired time of flight between
% Mars and Jupiter in days.
OPT.tof_total = [   300; ...
                    1500];      % tof w/ gravity assists [days]
                
% --- Set a margin on the time of flight windows between each planet
% In the first row, the first number is the lower margin and the second
% number is the upper margin. In this case, the user desires the trip to
% Mars to take between (300-50) and (300+100) days total.
OPT.tof_margin = [  50,100; ...
                    300,300];	% margin on tof -/+     [days]

% --- Launch Windows
% This is where the departure windows are set for the departure from
% "Earth". BOD.window1 sets the earilest launch date to be Jan 1, 2019 and
% BOD.window2 sets the latest launch date to be Jan 1, 2021.
BOD.window1 = [2019,1,1];       % earliest launch day	[year,month,day]
BOD.window2 = [2021,1,1];       % latest launch day  	[year,month,day]



%% Cost Function Selection

% --- Solver selection
% This is where the user selects which method is used to solve for the
% optimization problem. In this case, the Forward/Backward Shooting Method
% using an indirect solver is used. Another valid option for this
% optimization case is 'MGALT_DIR_FBSM_2D'.
% *The FSM will not work with multiple planet gravity assists, and STOpS
% will error out.
OPT.solver = 'MGALT_IN_FBSM_2D';



%% Optimization Options

% Note:
% To see the available options for all choices below, open the file titled 
% "Options" in the main directory, or right click on the function name 
% below and selecting "Open".


% --- Thrust Profile
% This is where the user can select which thrust profile to use for the
% spacecraft. Within the function, a user can add as many unique thrust 
% profiles as desired.
% For this demo, STOpS will use the same thrust profile as used in the Yam
% et al. verification cases.
OPT.thrust = optionsLowThrust(OPT,'Yam-STOUR');

% --- Island Parameters
% This is where the user can select the options for island connectivity,
% solution sharing, and number of migrations.
% For this demo, STOpS will use a case where four islands are used, each
% containing one of the optimization algorithms, and the islands will
% perform three migrations.
OPT.island = optionsIsland('3M-all');

% --- Cost and Weighting Parameters
% This is where the user can select the options for calculating the cost of
% each member. 
% For this demo, STOpS will use the default cost function values found by 
% Sheehan and Malloy to produce the best optimization results
[OPT.cost,OPT.weighting,OPT.ode] = optionsCost(OPT,'default');

% --- Algorithm Parameters
% This is where the user can select the options for each algorithm.
% *Note, if the choice for "optionsIsland" does not include one of the
%  algorithms listed below, it will not be used.
% For instance, if a secondary GA is desired, it would be called by
% defining the var as "OPT.GA(1,2) = parametersGA(...)".
% For this demo, STOpS will use a GA, DE, PSO, and MBH algorithm to explore
% the search space.
% GA = 100 population initial and 30 generations
% DE = 100 population initial and 30 generations
% PSO = 50 population initial and 75 time steps
% MBH = 3 migrations, 1200 population initial, and 500 iters on feasibility
OPT.GA(1,1) = parametersGA('100_30');           % Genetic Algorithm Island #1
OPT.DE(1,1) = parametersDE('100_30');           % Differential Evolution Island #1
OPT.PSO(1,1) = parametersPSO('50_75');          % Particle Swarm Island #1
OPT.MBH(1,1) = parametersMBH('3M-1200_500');	% Monotonic Basin Hopping Island #1



%% User Error Catch

% --- Check for errors
% The user does not need to make any choices here. This function will check
% all of the selections made above and if any of them will cause STOpS to
% error out later on, it will prevent the program from running right now.
% *To see this function in action, change "BOD.window2" to be [2025,1,1]
%  and STOpS will display an error message in the command window. 
if errorMain(BOD,OPT)
    return
end



%% Pre Optimization

% --- Start parallel processing
% If parallel processing is desired, this will start the parallel pool if
% it is not already initialized
if OPT.parallel
    if isempty(gcp())
        parpool('local')
    end
	OPT.ppool = gcp();
end

% --- Define the search window for the planetary vectors
% Calculate the earliest and latest Julian Dates for the whole search
% window and get the total time span for the optimization process
BOD.JD1 = getJulianDate(BOD.window1);           % earliest Julian date
BOD.JD2 = getJulianDate(BOD.window2) + sum(OPT.thrust.tt_end) + ...
    sum(OPT.thrust.time(:,end));             	% latest Julian Date;
BOD.tspan = [BOD.JD1,BOD.JD2];

% --- Define Position, Velocity, and Julian Date for Bodies
% Pulls necessary planetary data from "Data -> JPL_Horizons"
[BOD.bodies_R,BOD.bodies_V,BOD.bodies_JD] = planetBodies(BOD);

% --- Create variable limits
% Creates the limits used to govern the initial population generation. For
% instance, when the user selected a tof_margin of 50 and 100 for Earth to
% Mars with a desired tof of 300 days, this ensures that the tof variable 
% in the member array will be bounded between 250 and 400 days.
VAR = MGALT_varLimits(BOD,OPT);

% --- Get constants
% Pulls the necessary constants for the optimization process. By default,
% this function pulls all parameters for the sun, such as SOI, mu, etc...
% and other value constants such as AU and TU. This function will then work
% as a dynamic struct and will pull necessary constants for the planets
% selected.
% In this demo, this function will pull constants used to define Earth,
% Mars, and Jupiter.
CONST = getConstants(BOD);



%% Optimization

% --- Define selected as an empty struct
% Preallocation process
selected = struct;

% --- Timer
% Allow us to see how long the optimization process took out computer to
% perform
tic     % Start a timer

% --- Run through the solver
% This loop allows the solver to run for the desired number of migrations+1
for num_mig = 1:OPT.island.Nmig+1
    
    % --- Algorithm counter in case an algorithm is used more than once
    % For instance, if more than one GA was used, this number will index
    % the defined options between the first and second GA and allow STOpS
    % to known which options to use.
	count_GA = 1; 
    count_DE = 1; 
    count_PSO = 1; 
    count_MBH = 1;
    
    % --- Run through all the islands
    % This loop allows the solver to run through each island within the
    % island list
	for num_isl = 1:OPT.island.Nisl
        
        % Switch case to better organize the islands
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

% --- Get the total run time in minutes
run_time = toc/60;



%% Display Output

% ---  Display all info to the Command Window and plot trajectory
% This function will tell the user a lot of different things about the
% optimization process, including the actual tof between planets, the
% actual departure date, the amount of fuel used, etc... It will also plot
% the trajectory and the thruster pointing angles for the user to see
displayResults(BOD,CONST,OPT,VAR,run_time,eval_info)

% --- Open the file to load previous results
% Only for the demo script. Will open the file where a user can load
% already optimized trajectories
open Load_Results.m



%% Save File

% --- Get current path with pwd and set new path
% Get the user desired path from the options listed above
new_path = strcat(pwd,'\',OPT.save_folder);

% --- Get the solver name
% Shorten the name of the solver to make saving the file easier
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

% --- Temp name
% Generate a temp name by using today's date
counter = 1;
temp_name = strcat(date(),'_',solver_name,'_',num2str(counter),'.mat');

% --- Check to see if name already exists
% If the name already exists, then incriment the counter to make a new temp
% name
while isfile(fullfile(new_path,temp_name))
    counter = counter+1;
    temp_name = strcat(date(),'_',solver_name,'_',num2str(counter),'.mat');
end

% --- Save to the current path
% Will be in STOpS-MGALT
save(temp_name)

% --- Move results to the name of "save_folder"
% Move the saved data so it's not in the "STOpS-MGALT" folder and is
% instead in the desired save folder
movefile(temp_name,new_path);

% Message
fprintf('\n\n\nThe optimization process is complete!\n')
fprintf('The demo script has opened "Load_Results.m".\n')
fprintf('With "Load_Results", you can load previously run problems to generate plots and see the spacecraft performance.\n\n')
fprintf('Your file is named: \n"%s"\nand is located at: \n"%s"\n',temp_name,new_path);

% Clear vars
clear counter temp_name current_path new_path save_folder solver_name


