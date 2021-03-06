function [T] = SolarSail_thrust(OPT,sc_rad)
% FORM: [T] = SolarSail_thrust(OPT,sc_rad)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Thrust from a solar sail
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -sc_rad             (1,1)     	[float]         [AU]
% |         Radius in AU (scalar not vector)
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -T                  (1,1)       [float]         [N]
% |         Solar Sail Thrust
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -This Function was create to mimic the thrust function from 
% |     Near-Optimal Low-Thrust Orbit Transfers Generated by a 
% |     Genetic Algorithm
% |
% |-----------------------------------------------------------------------



%% Solve

% Constants
P0 = 49717.5705;    % Watts (kg*m^2/S^3) - Power avalable at 1 AU
g = 9.81;           % m/S^2

% Power
P = P0/sc_rad^2*(1.4279-(0.6139/sc_rad) + (0.0038/sc_rad^2))/(1 - 0.2619*sc_rad + 0.0797*sc_rad^2);
if P/P0 > 1.35
    P = 1.35*P0;
end

% Thrust
T = 2*P/g/OPT.thrust.Isp;   % Newtons



end


