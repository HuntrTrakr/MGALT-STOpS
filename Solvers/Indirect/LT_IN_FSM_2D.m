function [J,plot_vars] = LT_IN_FSM_2D(member,BOD,CONST,OPT,VAR)
% FORM: [J,plot_vars] = LT_IN_FSM_2D(member,BOD,CONST,OPT,VAR)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Cost function solver for use with MGALT STOpS
% |
% |     -This is the INDIRECT method for representing a low-thrust 
% |     orbital tragectory, taken from Conway
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -member             (1,Nvar)    [float]         [unitless]
% |         A single member of a population input into the function
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
% |     -CONST              (1,1)       [struct]        [unitless]
% |         A struct containing constants used in the calcs. Contains
% |         values for AU, TU, Sun (rad/mu/rp) and (rad/mu/rp/SOI/per) 
% |         for any bodies used in the optimization scheme. This is a 
% |         dynamic struct and will adapt to contain only the necesary 
% |         information
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -J                  (1,1)       [float]         [unitless]
% |     	The cost of this member, denoted as 'f' in other functions
% |     -plot_vars          (1,1)       [struct]     	[unitless]
% |         An object containing a lot of information about the 
% |         optimization parameters including: transfers(t and y ode 
% |         outputs), thrust values, thruster pointing angles, transfer 
% |         starting position, planet start/end locations for each 
% |         transfer, JD of each transfer, and tspans of each transfer
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -References
% |         B. A. Conway, Spacecraft Trajectory Optimization. Cambridge University Press, 2010.
% |
% |-----------------------------------------------------------------------



%% Setup

% Canonical Units
AU = CONST.AU;                  % [km/AU]
TU = CONST.TU*86400;            % [sec/TU]
mew_sun = CONST.Sun_mu;         % [km^3/s^2]
mew_sc = mew_sun*(TU^2/AU^3);  	% [DU^3/TU^2] Mew for the spacecraft

% Total Cost
J = 0;

% Other
transfers = VAR.transfers;



%% Pre-Allocation

% Orbit variables
tspan = zeros(transfers,2);
JD = zeros(transfers,2);


% Cost function variables
pos_rad_sc_trans = zeros(1,transfers);
pos_ang_sc_trans = zeros(1,transfers);
vel_rad_sc_trans = zeros(1,transfers);
vel_tan_sc_trans = zeros(1,transfers);
transfer_time_sc_trans = zeros(1,transfers);


% JD for the bodies during transfer segments
start = member(1,1);
tof = member(1,end);
JD(1,:) = [start, start+tof];                  % Start and end Julian Day

% tspan for all transfer segments
tspan(1,:) = [0, member(5)*86400]*(1/TU);	%TU



%% Initial Conditions               (Only for Departure Body)

% Get R and V for the departure planet at the first JD
[planet_R_dep,planet_V_dep,sc_dep_pos_rad,sc_dep_pos_ang,...
    sc_dep_vel_rad,sc_dep_vel_tan] = MGALT_conditionsInit(...
    JD(1,1),BOD,CONST,OPT,VAR,[1:3]);
sc_dep_mass = OPT.thrust.m0;      % inital mass (kg)

% For the plot variables
planet_departure = [planet_R_dep; planet_V_dep];

% Check to see if any additional dV from launch vehicle
if isfield(OPT.thrust,'launch_dV_rad')
    sc_dep_vel_rad = sc_dep_vel_rad + (OPT.thrust.launch_dV_rad*(TU/AU));
end
if isfield(OPT.thrust,'launch_dV_tan')
    sc_dep_vel_tan = sc_dep_vel_tan + (OPT.thrust.launch_dV_tan*(TU/AU));
end



%% Solve Intermediate Conditions    (Transfers Between)

% Solve between transfers
switch transfers
    
    case {1}        % If there are not any Gravity Assists
        
        % Define the state variable
        lam1 = member(2);
        lam2 = member(3);
        lam3 = member(4);
        Y0 = [lam1; lam2; lam3; sc_dep_pos_rad; sc_dep_pos_ang; sc_dep_vel_rad; sc_dep_vel_tan; sc_dep_mass];

        % ODE45 to solve for the orbit the with generated costates
        [ttot, Ytot] = ode45(@LT_IN_FSM_2D_EOM,...
            tspan,Y0,OPT.ode,CONST,OPT,mew_sc);

        % End Conditions
        pos_rad_sc_trans        = Ytot(end,4);
        pos_ang_sc_trans        = (Ytot(end,5)/360 - floor(Ytot(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_trans    	= Ytot(end,6);
        vel_tan_sc_trans     	= Ytot(end,7);
%         mass_sc_trans       	= Ytot(end,8);
        transfer_time_sc_trans  = ttot(end);

        % Append to spacecraft plotting variables
        plot_vars.Y0{1,1} = Y0;
        plot_vars.transfers{1,1} = [ttot,Ytot];
        
    otherwise       % If there are Gravity Assists
        
        errorPathDisplay()
        fprintf(2,'MGALT STOpS doe snot support FSM with more than two planets.\n')
        fprintf(2,'It is recommended to use the FBSM instead.\n\n')
        fprintf(2,'To see how Multiple FSM was implemented in STOpS, check versions before "MGALT Removal of Multiple FSM"\n')
        
end



%% Final Conditions                 (Only for Target Body)

% Get R and V for the target planet at the final JD
[planet_R_tar,planet_V_tar,sc_tar_pos_rad,sc_tar_pos_ang,sc_tar_vel_rad,...
    sc_tar_vel_tan] = MGALT_conditionsFinal(...
    JD(end,:),BOD,CONST,OPT,VAR,[4:6],pos_ang_sc_trans);

% For the plot variables
planet_target = [planet_R_tar; planet_V_tar];



%% Plotting Vars and Difference

% Plotting variables
plot_vars.planetary_conditions = [planet_departure, planet_target];
plot_vars.JD = JD;          % variables necessary to plot the thrust vectors
plot_vars.tspan = tspan;  	% time spans for the orbits

% plotOrbits(plot_vars,'MGALT_IN_FSM_2D',transfers)



%% Calculating Cost

% Get the total cost function for the population member
J = LT_FSM_costFun(J, ...
    [pos_rad_sc_trans, sc_tar_pos_rad], ...
    [pos_ang_sc_trans, sc_tar_pos_ang], ...
    [vel_rad_sc_trans, sc_tar_vel_rad], ...
    [vel_tan_sc_trans, sc_tar_vel_tan], ...
    transfer_time_sc_trans', ...
    CONST,OPT);



end


