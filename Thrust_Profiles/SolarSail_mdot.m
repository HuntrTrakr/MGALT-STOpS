function [mdot] = SolarSail_mdot(CONST,T,Isp)
% FORM: [mdot] = SolarSail_mdot(T,Isp)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -mdot from a solar sail
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -T                  (1,1)     	[float]         [N]
% |         Solar Sail Thrust
% |     -Isp                (1,1)      	[float]         [sec]
% |         Solar Sail Isp
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -mdot               (1,1)       [float]         [kg/s]
% |         Solar Sail mass flow rate
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -This Function was create to mimic the mdot function from 
% |     Near-Optimal Low-Thrust Orbit Transfers Generated by a 
% |     Genetic Algorithm
% |
% |-----------------------------------------------------------------------



%% Solve

% Constants
g = 9.81; % m/s^2

% Convert thrust
T = T/(CONST.TU*86400)^2*CONST.AU*1000;     % Newtons

% Mdot
mdot = T/g/Isp; % kg/s



end

