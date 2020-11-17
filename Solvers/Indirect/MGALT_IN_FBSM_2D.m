function [J,plot_vars] = MGALT_IN_FBSM_2D(member,BOD,CONST,OPT,VAR)
% FORM: [J,plot_vars] = MGALT_IN_FBSM_2D(member,BOD,CONST,OPT,VAR)
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

% ODE Parameters
% options = odeset('AbsTol',1e-6,'RelTol',1e-6);

% Other
transfers = VAR.transfers;



%% Pre-Allocation

% Orbit variables
tspan = zeros(transfers,2);
JD = zeros(transfers,2);


% ODE varaibles
tspan_divider = 400;


% Cost function variables
pos_rad_sc_fs = zeros(1,transfers);
pos_rad_sc_bs = zeros(1,transfers);
pos_ang_sc_fs = zeros(1,transfers);
pos_ang_sc_bs = zeros(1,transfers);
vel_rad_sc_fs = zeros(1,transfers);
vel_rad_sc_bs = zeros(1,transfers);
vel_tan_sc_fs = zeros(1,transfers);
vel_tan_sc_bs = zeros(1,transfers);
transfer_time_sc_fs = zeros(1,transfers);
transfer_time_sc_bs = zeros(1,transfers);


% Transfer variables
planet_departure = zeros(6,1);
planet_trans = zeros(6,transfers-1);
planet_target = zeros(6,1);



%% Perform Transfer(s)

switch transfers
    
    case {1}        % Going from planet A to B

        % JD for the bodies during transfer segments
        JD(1,:) = [member(1), member(1)+member(end)];     % Start and end Julian Day

        % tspan for all transfer segments
        tspan(1,:) = [0, member(end)*86400]*(1/TU);         % TU
        
        
        % ********** BODY CONDITIONS **********
        % Get departure body initial conditions
        [planet_R_dep,planet_V_dep,sc_dep_pos_rad,sc_dep_pos_ang,...
            sc_dep_vel_rad,sc_dep_vel_tan] = MGALT_conditionsInit(...
            JD(1,1),BOD,CONST,OPT,VAR,[1:3]);
        
        % Get target body final conditions
        [planet_R_tar,planet_V_tar,sc_tar_pos_rad,sc_tar_pos_ang,...
            sc_tar_vel_rad,sc_tar_vel_tan] = MGALT_conditionsInit(...
            JD(end,end),BOD,CONST,OPT,VAR,[4:6]);
        
        
        % ********** FORWARD SHOOTING **********
        % Check to see if any additional dV from launch vehicle
        if isfield(OPT.thrust,'launch_dV_rad')
            sc_dep_vel_rad = sc_dep_vel_rad + (OPT.thrust.launch_dV_rad*(TU/AU));
        end
        if isfield(OPT.thrust,'launch_dV_tan')
            sc_dep_vel_tan = sc_dep_vel_tan + (OPT.thrust.launch_dV_tan*(TU/AU));
        end
        
        % Define the coefficients for forward shooting
        lam1_fs = member(2);
        lam2_fs = member(3);
        lam3_fs = member(4);
        mass_fs_Y0 = OPT.thrust.m0;     % inital mass (kg)

        % Define the state varaible for forward shooting
        Y0_fs = [lam1_fs;lam2_fs;lam3_fs;sc_dep_pos_rad;sc_dep_pos_ang;...
            sc_dep_vel_rad;sc_dep_vel_tan;mass_fs_Y0];
        tspan_fs = linspace(tspan(1),tspan(end)/2,tspan_divider);

        % ODE45 to solve for forward segment
        [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc);
        
        % End Conditions
        pos_rad_sc_fs     	= Ytot_fs(end,4);
        pos_ang_sc_fs     	= (Ytot_fs(end,5)/360 - floor(Ytot_fs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_fs    	= Ytot_fs(end,6);
        vel_tan_sc_fs     	= Ytot_fs(end,7);
        transfer_time_sc_fs	= ttot_fs(end);
       
        
        % ********** BACKWARDS SHOOTING **********
        % Define the coefficients for backward shooting
        lam1_bs = member(5);
        lam2_bs = member(6);
        lam3_bs = member(7);
        
        % Assume a constant mass loss for the duration of the transfer
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs(1,end)-Ytot_fs(end,end));

        % Define the state varaible for backwards shooting
        Y0_bs = [lam1_bs;lam2_bs;lam3_bs;sc_tar_pos_rad;sc_tar_pos_ang;...
            sc_tar_vel_rad;sc_tar_vel_tan;mass_bs_Y0];
        tspan_bs = linspace(tspan(end)/2,tspan(1),tspan_divider);

        % ODE45 to solve for backwards segment
        [ttot_bs,Ytot_bs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_bs,Y0_bs,OPT.ode,CONST,OPT,mew_sc);
        
        % End Conditions
        pos_rad_sc_bs      	= Ytot_bs(end,4);
        pos_ang_sc_bs     	= (Ytot_bs(end,5)/360 - floor(Ytot_bs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_bs    	= Ytot_bs(end,6);
        vel_tan_sc_bs     	= Ytot_bs(end,7);
        transfer_time_sc_bs	= ttot_bs(1);

        
        % ********** Append to spacecraft plotting variables **********
        planet_departure = [planet_R_dep;planet_V_dep];
        planet_target = [planet_R_tar;planet_V_tar];
        plot_vars.Y0_fs{1,1} = Y0_fs;
        plot_vars.Y0_bs{1,1} = Y0_bs;
        plot_vars.transfers_fs{1,1} = [ttot_fs,Ytot_fs];
        plot_vars.transfers_bs{1,1} = [ttot_bs,Ytot_bs];
        plot_vars.transfers{1,1} = [ttot_fs,Ytot_fs;...
                                    flipud(ttot_bs)+ttot_fs(end),flipud(Ytot_bs)];
        
    otherwise       % Going from planet A to N via B, C, D, ...
        
        % Misc
        cost_count = 1;

        % For indexing the positions
        array_bodies = [4;5;6];
        array_member = [9,10];
        
        % JD for the bodies during transfer segments
        start = member(1);
        tof = member(11);
        JD(1,:) = [start, start+tof];
        for i1 = 2:transfers-1
            start = member(i1*11-10);
            tof = member(i1*11);
            JD(i1,:) = [start, start+tof];                  % Start and end Julian Day
        end
        start = member(end-7);
        tof = member(end);
        JD(end,:) = [start, start+tof];

        
        % tspan for all transfer segments
        tspan(1,:) = [0, member(11)*86400]*(1/TU);
        for i2 = 2:transfers-1
            tspan(i2,:) = [0, member(i2*11)*86400]*(1/TU);	%TU
        end
        tspan(end,:) = [0, member(end)*86400]*(1/TU);
        
        
        %% ********** DEPARTURE -> TRANSFER 1 **********
        
        % Departure body initial conditions
        [planet_R_dep,planet_V_dep,sc_dep_pos_rad,sc_dep_pos_ang,...
            sc_dep_vel_rad,sc_dep_vel_tan] = MGALT_conditionsInit(...
            JD(1,1),BOD,CONST,OPT,VAR,[1:3]);
        
        % Check to see if any additional dV from launch vehicle
        if isfield(OPT.thrust,'launch_dV_rad')
            sc_dep_vel_rad = sc_dep_vel_rad + (OPT.thrust.launch_dV_rad*(TU/AU));
        end
        if isfield(OPT.thrust,'launch_dV_tan')
            sc_dep_vel_tan = sc_dep_vel_tan + (OPT.thrust.launch_dV_tan*(TU/AU));
        end

        % Define the state varaible for forward shooting
        lam1_fs = member(2);
        lam2_fs = member(3);
        lam3_fs = member(4);
        mass_fs_Y0 = OPT.thrust.m0;     % inital mass (kg)
        
        Y0_fs = [lam1_fs;lam2_fs;lam3_fs;sc_dep_pos_rad;sc_dep_pos_ang;...
            sc_dep_vel_rad;sc_dep_vel_tan;mass_fs_Y0];
        tspan_fs = linspace(tspan(1,1),tspan(1,end)/2,tspan_divider);

        % ODE45 to solve for forward segment
        [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc);
        
        % FS end Conditions
        pos_rad_sc_fs(cost_count)     	= Ytot_fs(end,4);
        pos_ang_sc_fs(cost_count)     	= (Ytot_fs(end,5)/360 - floor(Ytot_fs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_fs(cost_count)    	= Ytot_fs(end,6);
        vel_tan_sc_fs(cost_count)     	= Ytot_fs(end,7);
        mass_sc_fs_end(cost_count)      = Ytot_fs(end,8);
        transfer_time_sc_fs(cost_count)	= ttot_fs(end);
        
        
        % Get transfer body position and s/c positions
        [planet_R_trans,planet_V_trans,sc_trans_pos_rad,sc_trans_pos_ang,...
            sc_trans_vel_rad,sc_trans_vel_tan,control,sc_vel_helio_enter] = ...
            MGALT_conditionsTransFBSM(...
            JD(1,2),...
            BOD,...
            CONST,...
            OPT,...
            VAR,...
            tspan(1,:),...
            mass_sc_fs_end,...
            [1,1],...
            [member(array_member),OPT.weighting.control_v],...
            sc_dep_pos_rad,...
            array_bodies);
        
        % Change the control numbers
        member(array_member) = control(1:2);

        % Define the state varaible for backwards shooting
        lam1_bs = member(5);
        lam2_bs = member(6);
        lam3_bs = member(7);
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs(1,end)-Ytot_fs(end,end));
        
        Y0_bs = [lam1_bs;lam2_bs;lam3_bs;sc_trans_pos_rad;sc_trans_pos_ang;...
            sc_trans_vel_rad;sc_trans_vel_tan;mass_bs_Y0];
        tspan_bs = linspace(tspan(1,end)/2,tspan(1,1),tspan_divider);

        % ODE45 to solve for backwards segment
        [ttot_bs,Ytot_bs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_bs,Y0_bs,OPT.ode,CONST,OPT,mew_sc);
        
        % BS end Conditions
        pos_rad_sc_bs(cost_count)      	= Ytot_bs(end,4);
        pos_ang_sc_bs(cost_count)     	= (Ytot_bs(end,5)/360 - floor(Ytot_bs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_bs(cost_count)    	= Ytot_bs(end,6);
        vel_tan_sc_bs(cost_count)     	= Ytot_bs(end,7);
        transfer_time_sc_bs(cost_count) = ttot_bs(1);
        

        % Count up the array bodies
        array_bodies = array_bodies+3;
        array_member = array_member+11;
        cost_count = cost_count+1;
        
        % Append to spacecraft plotting variables
        planet_departure = [planet_R_dep;planet_V_dep];
        planet_trans(:,1) = [planet_R_trans;planet_V_trans];    % Plot vars
        plot_vars.Y0_fs{1,1} = Y0_fs;
        plot_vars.Y0_bs{1,1} = Y0_bs;
        plot_vars.transfers_fs{1,1} = [ttot_fs,Ytot_fs];
        plot_vars.transfers_bs{1,1} = [ttot_bs,Ytot_bs];
        plot_vars.transfers{1,1} = [ttot_fs,Ytot_fs;...
                                    flipud(ttot_bs)+ttot_fs(end),flipud(Ytot_bs)];
        
        
        %%  ********** TRANSFER 1 -> TRANSFER n **********
        
        for i3 = 2:(transfers-1)
            
            % Calculate the FS departure position and vel
            rp_coe = member(array_member(1)-12);        % Number between 0 and 1 for flyby altitude
            
            % Get the gravity assist condition
            [sc_vel_helio_exit,~] = MGALT_conditionsGravityAssist(BOD,CONST,...
                planet_R_trans,planet_V_trans,sc_vel_helio_enter,rp_coe,i3);
            
            % Convert flyby LVLH
            [~,sc_dep_vel,~] = MGALT_convertHelioLVLH(...
                CONST,planet_R_trans,sc_vel_helio_exit);

            % Define the state varaible for forward shooting
            lam1_fs = member(end-6);
            lam2_fs = member(end-5);
            lam3_fs = member(end-4);
            mass_fs_Y0 = mass_bs_Y0;    % kg

            Y0_fs = [lam1_fs;lam2_fs;lam3_fs;sc_trans_pos_rad;sc_trans_pos_ang;...
                sc_dep_vel(1);sc_dep_vel(2);mass_fs_Y0];
            tspan_fs = linspace(tspan(i3,1),tspan(i3,end)/2,tspan_divider);

            % ODE45 to solve for forward segment
            [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
                tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc);

            % FS end Conditions
            pos_rad_sc_fs(cost_count)     	= Ytot_fs(end,4);
            pos_ang_sc_fs(cost_count)     	= (Ytot_fs(end,5)/360 - floor(Ytot_fs(end,5)/360))*360;   % Accounts for being larger than 360
            vel_rad_sc_fs(cost_count)    	= Ytot_fs(end,6);
            vel_tan_sc_fs(cost_count)     	= Ytot_fs(end,7);
            transfer_time_sc_fs(cost_count)	= ttot_fs(end);
            
            
            % Get transfer body position and s/c positions           
            [planet_R_trans,planet_V_trans,sc_trans_pos_rad,sc_trans_pos_ang,...
                sc_trans_vel_rad,sc_trans_vel_tan,control,sc_vel_helio_enter] = ...
                MGALT_conditionsTransFBSM(...
                JD(i3,2),...
                BOD,...
                CONST,...
                OPT,...
                VAR,...
                tspan(i3,:),...
                mass_sc_fs_end,...
                [1,1],...
                [member(array_member),OPT.weighting.control_v],...
                sc_dep_pos_rad,...
                array_bodies);

            % Change the control numbers
            member(array_member) = control(1:2);

            % Define the state varaible for backwards shooting
            lam1_bs = member(5);
            lam2_bs = member(6);
            lam3_bs = member(7);
            mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs(1,end)-Ytot_fs(end,end));

            Y0_bs = [lam1_bs;lam2_bs;lam3_bs;sc_trans_pos_rad;sc_trans_pos_ang;...
                sc_trans_vel_rad;sc_trans_vel_tan;mass_bs_Y0];
            tspan_bs = linspace(tspan(i3,end)/2,tspan(i3,1),tspan_divider);

            % ODE45 to solve for backwards segment
            [ttot_bs,Ytot_bs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
                tspan_bs,Y0_bs,OPT.ode,CONST,OPT,mew_sc);

            % BS end Conditions
            pos_rad_sc_bs(cost_count)      	= Ytot_bs(end,4);
            pos_ang_sc_bs(cost_count)     	= (Ytot_bs(end,5)/360 - floor(Ytot_bs(end,5)/360))*360;   % Accounts for being larger than 360
            vel_rad_sc_bs(cost_count)    	= Ytot_bs(end,6);
            vel_tan_sc_bs(cost_count)     	= Ytot_bs(end,7);
            transfer_time_sc_bs(cost_count) = ttot_bs(1);


            % Count up the array bodies
            array_bodies = array_bodies+3;
            array_member = array_member+11;
            cost_count = cost_count+1;

            % Append to spacecraft plotting variables
            planet_trans(:,i3) = [planet_R_trans;planet_V_trans];    % Plot vars
            plot_vars.Y0_fs{i3,1} = Y0_fs;
            plot_vars.Y0_bs{i3,1} = Y0_bs;
            plot_vars.transfers_fs{i3,1} = [ttot_fs,Ytot_fs];
            plot_vars.transfers_bs{i3,1} = [ttot_bs,Ytot_bs];
            plot_vars.transfers{i3,1} = [ttot_fs,Ytot_fs;...
                                        flipud(ttot_bs)+ttot_fs(end),flipud(Ytot_bs)];
            
        end
        
        
        %% ********** TRANSFER n -> TARGET **********
        
        % Calculate the FS departure position and vel
        rp_coe = member(array_member(1)-12);        % Number between 0 and 1 for flyby altitude
        
        % Get the gravity assist conditions
        [sc_vel_helio_exit,~] = MGALT_conditionsGravityAssist(BOD,CONST,...
            planet_R_trans,planet_V_trans,sc_vel_helio_enter,rp_coe,(size(BOD.bodies,2)-1));
        
        [~,sc_dep_vel,~] = MGALT_convertHelioLVLH(...
            CONST,planet_R_trans,sc_vel_helio_exit);
        
        % Define the state varaible for forward shooting
        lam1_fs = member(end-6);
        lam2_fs = member(end-5);
        lam3_fs = member(end-4);
        mass_fs_Y0 = mass_bs_Y0;    % kg
        
        Y0_fs = [lam1_fs;lam2_fs;lam3_fs;sc_trans_pos_rad;sc_trans_pos_ang;...
            sc_dep_vel(1);sc_dep_vel(2);mass_fs_Y0];
        tspan_fs = linspace(tspan(end,1),tspan(end,end)/2,tspan_divider);

        % ODE45 to solve for forward segment
        [ttot_fs,Ytot_fs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_fs,Y0_fs,OPT.ode,CONST,OPT,mew_sc);
        
        % FS end Conditions
        pos_rad_sc_fs(end)     	= Ytot_fs(end,4);
        pos_ang_sc_fs(end)     	= (Ytot_fs(end,5)/360 - floor(Ytot_fs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_fs(end)    	= Ytot_fs(end,6);
        vel_tan_sc_fs(end)     	= Ytot_fs(end,7);
        transfer_time_sc_fs(end)	= ttot_fs(end);
        
        
        % Get target body final conditions
        [planet_R_tar,planet_V_tar,sc_tar_pos_rad,sc_tar_pos_ang,...
            sc_tar_vel_rad,sc_tar_vel_tan] = MGALT_conditionsInit(...
            JD(end,end),BOD,CONST,OPT,VAR,array_bodies);
        
        
        % Define the coefficients for backward shooting
        lam1_bs = member(end-3);
        lam2_bs = member(end-2);
        lam3_bs = member(end-1);
        mass_bs_Y0 = mass_fs_Y0 - 2*(Ytot_fs(1,end)-Ytot_fs(end,end));

        % Define the state varaible for backwards shooting
        Y0_bs = [lam1_bs;lam2_bs;lam3_bs;sc_tar_pos_rad;sc_tar_pos_ang;...
            sc_tar_vel_rad;sc_tar_vel_tan;mass_bs_Y0];
        tspan_bs = linspace(tspan(end,end)/2,tspan(end,1),tspan_divider);

        % ODE45 to solve for backwards segment
        [ttot_bs,Ytot_bs] = ode45(@MGALT_IN_FBSM_2D_EOM,...
            tspan_bs,Y0_bs,OPT.ode,CONST,OPT,mew_sc);
        
        % End Conditions
        pos_rad_sc_bs(end)    	= Ytot_bs(end,4);
        pos_ang_sc_bs(end)     	= (Ytot_bs(end,5)/360 - floor(Ytot_bs(end,5)/360))*360;   % Accounts for being larger than 360
        vel_rad_sc_bs(end)    	= Ytot_bs(end,6);
        vel_tan_sc_bs(end)     	= Ytot_bs(end,7);
        transfer_time_sc_bs(end)	= ttot_bs(1);
        
        % Plot Vars
        planet_target = [planet_R_tar;planet_V_tar];
        plot_vars.Y0_fs{end+1,1} = Y0_fs;
        plot_vars.Y0_bs{end+1,1} = Y0_bs;
        plot_vars.transfers_fs{end+1,1} = [ttot_fs,Ytot_fs];
        plot_vars.transfers_bs{end+1,1} = [ttot_bs,Ytot_bs];
        plot_vars.transfers{end+1,1} = [ttot_fs,Ytot_fs;...
                                    flipud(ttot_bs)+ttot_fs(end),flipud(Ytot_bs)];

                                
end



%% Plotting Vars and Difference

% Plotting variables
plot_vars.planetary_conditions = [planet_departure, planet_trans, planet_target];
plot_vars.JD = JD;          % variables necessary to plot the thrust vectors
plot_vars.tspan = tspan;  	% time spans for the orbits

% plotOrbits(plot_vars,'MGALT_IN_FBSM_2D',transfers)



%% Calculating Cost

% Get the total cost function for the population member
J = MGALT_FBSM_costFun(J, ...
    [pos_rad_sc_fs', pos_rad_sc_bs'], ...
    [pos_ang_sc_fs', pos_ang_sc_bs'], ...
    [vel_rad_sc_fs', vel_rad_sc_bs'], ...
    [vel_tan_sc_fs', vel_tan_sc_bs'], ...
    [transfer_time_sc_fs', transfer_time_sc_bs'],...
    CONST,OPT,VAR);



end


