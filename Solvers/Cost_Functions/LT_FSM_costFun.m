function [J_fin] = LT_FSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,...
    CONST,OPT)
% FORM: [J_fin] = LT_FSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,...
%       CONST,OPT)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function calculates the cost of a FSM transfer between two 
% |     planets and takes into account the different selection parameters 
% |     for the cost funciton
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -J_init             (1,1)       [float]         [unitless]
% |     	The cost of this member, denoted as 'f' in other functions
% |     -pos        	(2,transfers) 	[float]      	[DU]
% |         The radial positions of the spacecraft and planet (sc;planet)
% |     -ang        	(2,transfers) 	[float]      	[deg]
% |         The angulat positions of the spacecraft and planet (sc;planet)
% |     -vel_rad       	(2,transfers) 	[float]      	[DU/TU]
% |         The radial velocity of the spacecraft and planet (sc;planet)
% |     -vel_tan       	(2,transfers) 	[float]      	[DU/TU]
% |         The angulat velocity of the spacecraft and planet (sc;planet)
% |     -tof        	(2,transfers) 	[float]      	[TU]
% |         The times of flight for each transfer
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
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -J_fin              (1,1)       [float]         [unitless]
% |     	The cost of this member, denoted as 'f' in other functions
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initials

% Constants
J_cost = 0;



%% Cost
   
% ---------- Minimize Final Position Error ----------------------------
JR = (pos(1,1) - pos(1,2))^2/(OPT.cost.tolR^2);
J_cost = J_cost + (OPT.cost.R*JR);

% ---------- Minimize Final Angular Displacement Error ----------------
JT = (ang(1,1) - ang(1,2))^2/(OPT.cost.tolTheta^2);
J_cost = J_cost + (OPT.cost.Theta*JT);

% ---------- Minimize Final Radial Velocity Error ---------------------
JU = (vel_rad(1,1) - vel_rad(1,2))^2/(OPT.cost.tolU^2);
J_cost = J_cost + (OPT.cost.U*JU);

% ---------- Minimize Final Tangential Velocity Error -----------------
JV = (vel_tan(1,1) - vel_tan(1,2))^2/(OPT.cost.tolV^2);
J_cost = J_cost + (OPT.cost.V*JV);

% ---------- Restraint on Time of Flight ------------------------------
TU = CONST.TU*86400;                     % [sec/TU]
tt_end = OPT.thrust.tt_end(1)*86400*(1/TU);    % TU
Jtt = (tof(1,1) - tt_end)^2/(OPT.weighting.W_tof_conv^2);
J_cost = J_cost + (OPT.cost.tt*Jtt);


% Add all the costs together
J_fin = J_cost + J_init;



end


