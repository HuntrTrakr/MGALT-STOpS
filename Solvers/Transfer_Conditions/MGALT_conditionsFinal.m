function [R_fin_return,V_fin_return,pos_rad,pos_ang,vel_rad,vel_tan] = ...
    MGALT_conditionsFinal(JD,BOD,CONST,OPT,VAR,planet_parse,arc_trans)
% FORM: [R_fin_return,V_fin_return,pos_rad,pos_ang,vel_rad,vel_tan] = ...
%       MGALT_conditionsFinal(JD,bodies_R,bodies_V,bodies_JD,theta,arc_trans)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function gets a lot of different parameters for the target 
% |     planet. These include the R and V of the target at the final 
% |     Julian Date, as well as the radial position and velocity and the 
% |     angular location and velocity
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -JD                 (1,1)       [float]         [JD]
% |         Julian day for calculation
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
% |     -planet_parse       (3,1)       [int]        [unitless]
% |         An array containing the index numbers for the bodies
% |     -arc_trans        	(1,1)   	[float]         [deg]
% |         The planet angular position
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -R_fin_return       (3,1)       [float]         [AU]
% |         Radius vector
% |     -V_fin_return       (3,1)       [float]         [AU/TU]
% |         Velocity vector
% |     -pos_rad        	(1,1)   	[float]         [DU]
% |         The planet radial position
% |     -pos_ang        	(1,1)   	[float]         [deg]
% |         The planet angular position
% |     -vel_rad        	(1,1)   	[float]         [DU/TU]
% |         The planet radial velocity
% |     -vel_tan        	(1,1)   	[float]         [DU/TU]
% |         The planet tangential velocity
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Solve final spacecraft conditions

% Initials
AU = CONST.AU;         	% [km/AU]
TU = CONST.TU*86400;	% [sec/TU]


% For planet intercept
if OPT.cost.Theta
    
    % Get R and V for the transfer planet to the target planet
    [R_fin_return,V_fin_return] = MGALT_stateBodies(JD(2),BOD,CONST,OPT,VAR,planet_parse);
    
    % Transform the R and V by the distance units
    R_fin = R_fin_return*(1/AU);    % DU
    V_fin = V_fin_return*(TU/AU);   % TU
    
    % FUTURE WORK
    %{
    Removing the ability for heliocentric revolutions, as the transfer
    which these revolutions take place on would need to be specified and
    that's not within the scope of my work at the moment. If time permits,
    I'll come back and try to impliment this into the cost functions as a
    random generation into the initial population member.
    
    This code snippit was taken from the "EP_cost_fun_segmented_2D.m" file
    on line 156 from Sheehan's version of STOpS.
    
    "
    if strcmp(thrust_opt.orbit_check,'on')
        orbits = thrust_opt.orbits; 
        arc_des = arcf_test - arc0;
        if arc_des <= 0
            arc_des = arc_des + 360;
        end
        arc_des = arc_des + 360*(orbits-1);
        arcf_test = arc_des + arc0;
    end
    "
    
    %}


    
% For orbit intercept
else
    
    % Get the position where the planet is at the first JD
    [R_planet_start,V_planet_start] = MGALT_stateBodies(JD(1),BOD,CONST,OPT,VAR,planet_parse);
    Y0_planet_start = [R_planet_start;V_planet_start];
    
    % Get the period for the planet's whole orbit
    planet_period = CONST.(strcat(BOD.bodies{end},"_per"));     % https://blogs.mathworks.com/loren/2005/12/13/use-dynamic-field-references/
    planet_tspan = [0, planet_period*86400];     %s
    
    % Get the location of the planet over the entirity of its orbit
    [~,planet_data] = ode45(@orbit3D,planet_tspan,Y0_planet_start,OPT.ode,CONST.Sun_mu);
    
    % Get the R and V positions with the angle
    planet_R_orbit = planet_data(:,1:3);
    planet_V_orbit = planet_data(:,4:6);
    planet_angle_orbit = atan2d(planet_R_orbit(:,2),planet_R_orbit(:,1)); % (deg) ccw from positive X
    for i1 = 1:length(planet_angle_orbit)
        if planet_angle_orbit(i1) < 0
            planet_angle_orbit(i1) = 360 + planet_angle_orbit(i1);
        end
    end
    
    % Find where the minimum is located on the transfer
    [~, index] = min(abs(planet_angle_orbit-arc_trans));
    R_orbit = planet_R_orbit(index,:)'; % DU
    V_orbit = planet_V_orbit(index,:)'; % TU


    % Transform the R and V by the distance units
    R_fin = R_orbit*(1/AU);    % DU
    V_fin = V_orbit*(TU/AU);   % TU
    
    % Get R and V for the target planet at the final position
    [R_fin_return,V_fin_return] = MGALT_stateBodies(JD(2),BOD,CONST,OPT,VAR,planet_parse);
    
end
    
% Convert into polar coords
pos_rad = norm(R_fin);                  % DU
pos_ang = atan2d(R_fin(2),R_fin(1));	% (deg) ccw from positive X

% Account for tangent region
if pos_ang < 0
    pos_ang = 360 + pos_ang;
end

% Get the coords in LVLH
Rf_lvlh = R_fin/pos_rad;
Wf_lvlh = cross(R_fin,V_fin)/norm(cross(R_fin,V_fin));
Sf_lvlh = cross(Wf_lvlh,Rf_lvlh);
Tf_ECI_lvlh = [Rf_lvlh, Sf_lvlh, Wf_lvlh];
Vf_lvlh = Tf_ECI_lvlh'*V_fin;      % DU/TU

% Radial and Tangential Velocity for the solver
vel_rad = Vf_lvlh(1);           % DU/TU, parallel to R vector
vel_tan = Vf_lvlh(2);           % DU/TU, perpendicular to R vector



end


