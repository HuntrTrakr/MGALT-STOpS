function [vel_sc_exit,rp] = MGALT_conditionsGravityAssist(BOD,CONST,...
    pos_body,vel_body,vel_sc_enter,coe,index)
% FORM: [vel_sc_exit,rp] = MGALT_conditionsGravityAssist(BOD,CONST,...
%       pos_body,vel_body,vel_sc_enter,coe,index)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function gets the spacecraft velocity after a gravity assist,
% |     as well as the radius of perigee
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
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
% |     -pos_body           (3,1)    	[float]         [AU]
% |         Heliocentric position of the planet
% |     -vel_body           (3,1)       [float]         [AU/TU]
% |         Heliocentric velocity of the planet
% |     -vel_sc_enter      	(3,1)    	[float]         [AU/TU]
% |         Heliocentric velocity of the spacecraft before gravity assist
% |     -coe                (1,1)       [float]         [unitless]
% |         Coefficient to determine rp value
% |     -index            	(1,1)       [int]           [unitless]
% |         Used to determine current planet
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -vel_sc_exit      	(3,1)       [float]         [DU/TU]
% |         Heliocentric radial velocity of spacecraft after flyby
% |     -rp                 (1,1)       [float]         [km]
% |         Radius of perigee for spacecraft flyby
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initial

% Some basic info
pos_center = [0;0;0];
mue_sun = CONST.Sun_mu;                                 % km^3/s^2
mue = CONST.(strcat(BOD.bodies{index},"_mu"));       	% km^3/s^2
SOI = CONST.(strcat(BOD.bodies{index},"_SOI"))/2;       % km *radius
min_flyby = CONST.(strcat(BOD.bodies{index},"_rp"));	% km *radius



%% Solve for the Exit Velocity Vectors

% Solve for the unit vectors, Curtis: 8.75
u_v = vel_body/norm(vel_body);              % planet's velocity vector
planet_to_sun = pos_center - pos_body;
u_s = planet_to_sun/norm(planet_to_sun);	% planet pointing to the sun


% Get the scalar components of V1_vec, Curtis: 8.76
% Angle between V1_vec and V_vec, Curtis: 8.76
cos_alpha = dot(vel_sc_enter,vel_body)/(norm(vel_sc_enter)*norm(vel_body));
alpha_1 = acos(cos_alpha);                  % rad


% Calculate V, Curtis: 8.79
V = sqrt(mue_sun/norm(pos_body));


% Calculate V_inf1 vec, Curtis: 8.72
V_inf1_vec = vel_sc_enter-vel_body;
V_inf1_mag = norm(V_inf1_vec);


% Get the scalar components of V_inf1_vec, Curtis: 8.81
V_inf1_V = norm(vel_sc_enter)*cos(alpha_1) - V;
V_inf1_S = norm(vel_sc_enter)*sin(alpha_1);


% Angle between V_inf1_vec and V_vec, Curtis: 8.84
phi_1 = atan2(V_inf1_S,V_inf1_V);   % rad
if phi_1 < 0    % Quadrant ambiguity
    phi_1 = (2*pi)+phi_1;           % rad
end


% Calculate the turning angle, delta, Curtis: 8.54
% Get v_inf_mag, Curtis: 8.50 and 8.51
if norm(vel_sc_enter) <= norm(vel_body)
    v_inf_mag = norm(vel_body)-norm(vel_sc_enter);
else
    v_inf_mag = norm(vel_sc_enter)-norm(vel_body);
end
rp = coe*(SOI-min_flyby)+min_flyby;  % Get radius perigee from un-normed data 
delta = 2*asin(1/(1+((rp*(v_inf_mag^2))/mue)));     % rad


% Angle between V_inf2_vec and V_vec, Curtis: 8.85
phi_2 = phi_1+delta;	% rad


% Calculate V_inf2_vec, Curtis: 8.86
V_inf2_vec = V_inf1_mag*cos(phi_2)*u_v + V_inf1_mag*sin(phi_2)*u_s;


% V2 vec, Curtis: 8.87
vel_sc_exit = vel_body+V_inf2_vec;



end


