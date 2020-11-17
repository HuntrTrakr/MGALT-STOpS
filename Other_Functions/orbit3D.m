function [dY] = orbit3D(t,Y,mew)
% FORM: [dY] = orbit3D(t,Y,mew)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Return R, V, and JD arrays for bodies included in the problem for 
% |     the valid time span
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -t                  (1,1)     	[float]       	[s]
% |         Cell containing strings of the planets
% |     -Y                  (6,1)     	[float]         [km][km/s]
% |         State vector, see MISC for unpacking
% |     -mew                (1,1)     	[float]         [km^3/s^2]
% |         Gravitaional parameter of the sun
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -dY                 (6,1)       [float]         [km][km/s]
% |         Derivative of the state vector
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -Unpacking
% |         Y(1) = Rx       [km]
% |         Y(2) = Ry      	[km]
% |         Y(3) = Rz   	[km]
% |         Y(4) = Vx       [km/s]
% |         Y(5) = Vy       [km/s]
% |         Y(6) = Vz       [km/s]
% |
% |     -Used:
% |         [time,dY] = ode45(@orbit3D,tspan,Y,options,mew)
% |
% |-----------------------------------------------------------------------



%% EOM's

% Unpack
Rx = Y(1);          % [km]
Ry = Y(2);          % [km]
Rz = Y(3);          % [km]
R = [Rx; Ry; Rz]; 	% [km]
r = norm(R);        % [km]
Vx = Y(4);          % [km/s]
Vy = Y(5);          % [km/s]
Vz = Y(6);          % [km/s]

% Derivatives
dR = [Vx; Vy; Vz];  % [km/s]
dV = -mew*R/r^3;    % [km^3/s^2]

% Output
dY = [dR;dV];



end


