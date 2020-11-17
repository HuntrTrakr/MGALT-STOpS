function [sc_R0_lvlh,sc_V0_lvlh,sc_T0_ECI_lvlh] = ...
    MGALT_convertHelioLVLH(CONST,sc_R0_helio,sc_V0_helio)
% FORM: [sc_R0_lvlh,sc_V0_lvlh,sc_T0_ECI_lvlh] = ...
%       MGALT_convertHelioLVLH(CONST,sc_R0_helio,sc_V0_helio)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function converts a spacecraft's initial R and V vectors into
% |     heliocentric rad and tan vectors
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -sc_R_init          (3,1)       [float]         [AU]
% |         Radius vector of the current planet
% |     -sc_V_init          (3,1)    	[float]         [AU]
% |         Heliocentric velocity of the spacecraft after gravity assist
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -sc_vel_rad         (3,1)       [float]         [DU/TU]
% |         Heliocentric radial velocity component
% |     -sc_vel_tan         (3,1)       [float]         [DU/TU]
% |         Heliocentric tangential velocity component
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Convert

% Constants
AU = CONST.AU;                  % [km/AU]
TU = CONST.TU*86400;            % [sec/TU]

% R and V in AU and TU
sc_R0 = sc_R0_helio * (1/AU); 	% [AU]
sc_V0 = sc_V0_helio * (TU/AU);	% [TU]

% Convert into polar coords from cartesian
sc_pos_rad = norm(sc_R0);               % [AU]

% Get the coords in LVLH
sc_R0_lvlh = sc_R0/sc_pos_rad;          % [AU]
sc_W0_lvlh = cross(sc_R0,sc_V0)/norm(cross(sc_R0,sc_V0));
sc_S0_lvlh = cross(sc_W0_lvlh,sc_R0_lvlh);
sc_T0_ECI_lvlh = [sc_R0_lvlh, sc_S0_lvlh, sc_W0_lvlh];
sc_V0_lvlh = sc_T0_ECI_lvlh'*sc_V0;     % [AU/TU]



end


