function [planet_R_init,planet_V_init,...
    sc_pos_rad,sc_pos_ang,sc_vel_rad,sc_vel_tan,control,sc_V_init] = ...
    MGALT_conditionsTransFBSM(JD,BOD,CONST,OPT,VAR,tspan,mass,...
    thrust_segments,control,pos_rad,array_bodies)
% FORM: [planet_R_init,planet_V_init,...
%       sc_pos_rad,sc_pos_ang,sc_vel_rad,sc_vel_tan,control,sc_V_init] = ...
%       MGALT_conditionsTransFBSM(JD,BOD,CONST,OPT,VAR,tspan,mass,...
%       thrust_segments,control,pos_rad,array_bodies)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function gets a lot of different parameters for the 
% |     departure planet. These include the R and V of the departure at 
% |     the first Julian Date, as well as the radial position and 
% |     velocity and the angular location and velocity
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -JD                 (1,1)       [float]         [JD]
% |         Julian day for calculation at departure
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
% |     -mass               (1,1)       [float]         [kg]
% |         Mass of the spacecraft upon arrival
% |     -thrust_segments	(1,n)       [float]         [AU]
% |         The thrusting segments for the direct and indirect method
% |     -control            (1,3)       [float]         [unitless]
% |         The three control variables. The X/Y position are (1)/(2) and
% |         the tolerance is (3)
% |     -pos_rad            (1,1)       [float]         [AU]
% |         The radial position of the S/C with respect to the sun
% |     -array_bodies       (3,1)       [int]        [unitless]
% |         An array containing the index numbers for the bodies
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -planet_R_init    	(3,1)       [float]         [AU]
% |         Radius vector for the planet
% |     -planet_V_init     	(3,1)       [float]         [AU/TU]
% |         Velocity vector for the planet
% |     -sc_pos_rad        	(1,1)   	[float]         [DU]
% |         The planet radial position
% |     -sc_pos_ang        	(1,1)   	[float]         [deg]
% |         The planet angular position
% |     -sc_vel_rad        	(1,1)   	[float]         [DU/TU]
% |         The planet radial velocity
% |     -sc_vel_tan        	(1,1)   	[float]         [DU/TU]
% |         The planet tangential velocity
% |     -control            (1,3)       [float]         [unitless]
% |         The three control variables. The X/Y position are (1)/(2) and
% |         the tolerance is (3)
% |     -sc_V_init        	(3,1)   	[float]         [DU/TU]
% |         Velocity vector for the spacecraft
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initials

% Orbital stuff
AU = CONST.AU;        	% [km/AU]
TU = CONST.TU*86400;	% [sec/TU]

% S/C control stuff
control_x = control(1);
control_y = control(2);
tol = control(3);
t_init = tspan(1)*TU;   % sec
t_final = tspan(2)*TU;  % sec

% S/C info
n_available = OPT.thrust.n_available;
time_step = 1;      % Assuming 1 because not a true SFT where each time step has its own constraints, only stepping between planets



%% Planet

% Get R and V for the specified planet at the specified JD
[planet_R_init,planet_V_init] = MGALT_stateBodies(JD,BOD,CONST,OPT,VAR,array_bodies);



%% Thrust and Duty Cycle

% Get the duty cycle bsed off the type
switch OPT.thrust.duty_cycle_type
    
    case {'calculated'}
        
        % Look at the different thrust methods
        switch OPT.thrust.thrust_method
            
            case {'constant'}
                
                duty_cycle = sum(thrust_segments)/size(thrust_segments,2);
                
            case {'variable'}
                
                duty_cycle = 1; % Assume always on
                
            case {'equation'}
                
                % Assume on or off, use ODE later
                duty_cycle = sum(thrust_segments)/size(thrust_segments,2);
                
            otherwise
                
                errorPathDisplay();
                fprintf('Incorrect thrust profile selected.\n\n')
                return
                
        end
        
    case {'constant'}
        
        % Already defined as a number in OPT.thrust.duty_cycle, acts as a
        % pass
        duty_cycle = OPT.thrust.duty_cycle;
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect choice for the duty cycle type.\n')
        return
        
end

% Get the thrust based off the type
switch OPT.thrust.thrust_method
    
    case {'constant'}
        
        thrust = OPT.thrust.thrust;
        
    case {'variable'}
        
        thrust = max(OPT.thrust.thrust);
        
    case {'equation'}
        
        % Incorrect assumption to use the thrust at departure planet
        thrust = feval(OPT.thrust.thrust,OPT,pos_rad);	% N
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect choice for the thrust type.\n')
        return
        
end

thrust_max = thrust/n_available;



%% S/C

sc_R_init = planet_R_init;                  % Assume the s/c is leaving from the planet center

% Calculate the max velocity for the s/c	% Englander, Automated Solutions (2)
dV_max = ((duty_cycle*n_available*thrust_max)*(t_final-t_init))/(mass*time_step);      % m/s
dV_max = dV_max/1000;                       % km/s

% Calculate the change in velocity from the s/c
vel_norm = norm([control_x,control_y]);
if vel_norm > 1
    
    while vel_norm > 1
        if control_x > 0
            control_x = control_x - tol;
        else
                control_x = control_x + tol;
        end
        if control_y > 0
            control_y = control_y - tol;
        else
                control_y = control_y + tol;
        end
        vel_norm = norm([control_x,control_y]);
    end
    
end
add_x_dV = control_x*dV_max;
add_y_dV = control_y*dV_max;

% Add into current planetary position
sc_V_init = planet_V_init;
sc_V_init(1) = planet_V_init(1)+add_x_dV;
sc_V_init(2) = planet_V_init(2)+add_y_dV;


% R and V in DU and TU
sc_R0 = sc_R_init*(1/AU);             % DU
sc_V0 = sc_V_init*(TU/AU);            % TU

% Convert into polar coords
sc_pos_rad = norm(sc_R0);             % DU
sc_pos_ang = atan2d(sc_R0(2),sc_R0(1));	% (deg) ccw from positive X

% Account for tangent region
if sc_pos_ang < 0
    sc_pos_ang = 360 + sc_pos_ang;
end

% Get the coords in LVLH
sc_R0_lvlh = sc_R0/sc_pos_rad;
sc_W0_lvlh = cross(sc_R0,sc_V0)/norm(cross(sc_R0,sc_V0));
sc_S0_lvlh = cross(sc_W0_lvlh,sc_R0_lvlh);
sc_T0_ECI_lvlh = [sc_R0_lvlh, sc_S0_lvlh, sc_W0_lvlh];
sc_V0_lvlh = sc_T0_ECI_lvlh'*sc_V0;      % DU/TU

% Radial and Tangential Velocity for the solver
sc_vel_rad = sc_V0_lvlh(1);           % DU/TU, parallel to R vector
sc_vel_tan = sc_V0_lvlh(2);           % DU/TU, perpendicular to R vector

% Control
control(1) = control_x;
control(2) = control_y;




end


