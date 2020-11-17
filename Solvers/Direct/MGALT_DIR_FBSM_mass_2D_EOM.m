function [d_mass] = MGALT_DIR_FBSM_mass_2D_EOM(t,mass,CONST,OPT,pos_rad,mem_thrust)
% FORM: [dY] = MGALT_DIR_FBSM_2D_EOM(t,Y,AU,TU,mew,thrust,phi)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Equations of motion to get the change in mass over a certain time
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -t                  (1,1)       [float]         [unitless]
% |         Time, unused
% |     -mass               (1,1)       [float]         [kg]
% |         Mass of the spacecraft initial
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
% |     -pos_rad            (1,1)       [float]         [AU]
% |         The radial position of the S/C with respect to the sun
% |     -mem_thrust         (1,Nseg)	[bool][float] 	[unitless][N]
% |         Binary indicator if the s/c is thrusting or not
% |         The thrust for each angle
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -d_mass             (1,1)       [float]         [kg]
% |         change in the spacecraft's mass
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Get Mass

% Interpret if SC is thrusting and what the thrust is
T = getSCThrust(CONST,OPT,pos_rad,mem_thrust);   % (kg*DU/TU^2)

% Derivatives
d_mass = -1*getSCmdot(CONST,OPT,T);       % kg/TU



end


