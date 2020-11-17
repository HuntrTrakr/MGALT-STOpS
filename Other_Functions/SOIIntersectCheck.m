function [does_intersect,body_pos,values] = ...
    SOIIntersectCheck(BOD,CONST,OPT,VAR,plot_vars,index)
% FORM: [does_intersect,body_pos,values] = ...
%       SOIIntersectCheck(BOD,CONST,OPT,VAR,plot_vars,index)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function acts as a wrapper for "SOIIntersectLocation" 
% |     by extracting the planet string, finding the body_pos, and 
% |     providing it with the s/c X,Y,vX,vY components by parsing out 
% |     the plot_vars
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
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -plot_vars          (1,1)       [struct]     	[unitless]
% |         An object containing a lot of information about the 
% |         optimization parameters including: transfers(t and y ode 
% |         outputs), thrust values, thruster pointing angles, transfer 
% |         starting position, planet start/end locations for each 
% |         transfer, JD of each transfer, and tspans of each transfer
% |     -index              (1,1)       [int]       	[unitless]
% |         Index for the current body number
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -does_intersect 	(1,1)       [boolean]    	[unitless]
% |         If the s/c trajectory intersects the target SOI at any point
% |         along the orbital path
% |     -body_pos        	(1,1)       [float]         [AU]
% |         Heliocentric X,Y components of the desired transfer body at a
% |         predetermined JD
% |     -values             (2,2)       [float]         [AU][AU/TU]
% |         Position (X,Y); Velocity (X,Y) for the spacecraft
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Solve

% Find where the target body is at the target JD
[body_pos,~] = MGALT_stateBodies(plot_vars.JD(1,end),BOD,CONST,OPT,VAR,[4:6]);
% *Hardcoded as 4:6 because only used in the FSM which is only 2 planets


% Different methods for solver
switch OPT.solver

    case {'LT_DIR_FSM_2D'}

        sc_pos_rad = plot_vars.transfers{1,1}(:,2);     	% (DU) radial position
        sc_pos_ang = plot_vars.transfers{1,1}(:,3).*pi/180;	% (rad) angular position
        sc_vel_rad = plot_vars.transfers{1,1}(:,4);         % (DU/TU) radial velocity
        sc_vel_ang = plot_vars.transfers{1,1}(:,5);         % (DU/TU) tangential velocity
        
    case {'LT_IN_FSM_2D'}
        
        sc_pos_rad = plot_vars.transfers{1,1}(:,5);         % (DU) radial position
        sc_pos_ang = plot_vars.transfers{1,1}(:,6)*pi/180;	% (rad) angular position
        sc_vel_rad = plot_vars.transfers{1,1}(:,7);         % (DU/TU) radial velocity
        sc_vel_ang = plot_vars.transfers{1,1}(:,8);         % (DU/TU) tangential velocity

    otherwise

        errorPathDisplay();
        errorSolver();
        return

end

[sc_X,sc_Y] = pol2cart(sc_pos_ang,sc_pos_rad);          % Both in DU
[sc_vX,sc_vY] = pol2cart(sc_vel_ang,sc_vel_rad);        % Both in DU
        
[does_intersect,values] = SOIIntersectLocation(BOD,CONST,body_pos,[sc_X,sc_Y],[sc_vX,sc_vY],index);



end


