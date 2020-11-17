function [J_fin] = MGALT_FBSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,...
    CONST,OPT,VAR)
% FORM: [J_fin] = MGALT_FBSM_costFun(J_init,pos,ang,vel_rad,vel_tan,tof,...
%       CONST,OPT,VAR)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function calculates the cost of a particular member by 
% |     looking at the results from the planetary transfers and comparing 
% |     them to the desired conditions
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
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
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
array_bodies = [4;5;6];
J_trans = zeros(1,VAR.transfers);



%% Cost

for i4 = 1:VAR.transfers
   
    % ---------- Minimize Final Position Error ----------------------------
    JR = (pos(i4,1) - pos(i4,2))^2/(OPT.cost.tolR^2);
    J_trans(i4) = J_trans(i4) + (OPT.cost.R*JR);

    % ---------- Minimize Final Angular Displacement Error ----------------
    JT = (ang(i4,1) - ang(i4,2))^2/(OPT.cost.tolTheta^2);
    J_trans(i4) = J_trans(i4) + (OPT.cost.Theta*JT);

    % ---------- Minimize Final Radial Velocity Error ---------------------
    JU = (vel_rad(i4,1) - vel_rad(i4,2))^2/(OPT.cost.tolU^2);
    J_trans(i4) = J_trans(i4) + (OPT.cost.U*JU);

    % ---------- Minimize Final Tangential Velocity Error -----------------
    JV = (vel_tan(i4,1) - vel_tan(i4,2))^2/(OPT.cost.tolV^2);
    J_trans(i4) = J_trans(i4) + (OPT.cost.V*JV);

    % ---------- Restraint on Time of Flight ------------------------------
    TU = CONST.TU*86400;                     % [sec/TU]
    tt_end = OPT.thrust.tt_end(i4)*86400*(1/TU);    % TU
    Jtt = (tof(i4,1) + tof(i4,2) - tt_end)^2/(OPT.weighting.W_tof_conv^2);
    J_trans(i4) = J_trans(i4) + (OPT.cost.tt*Jtt);
    
    array_bodies = array_bodies+3;
    
end

% Add all the costs together
J_fin = sum(J_trans) + J_init;



end


