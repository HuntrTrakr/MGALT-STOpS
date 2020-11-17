function [sc_R0_helio,sc_V0_helio] = ...
    MGALT_convertLVLHHelio(CONST,sc_R0_lvlh,sc_V0_lvlh,sc_T0_ECI_lvlh)
% FORM: [sc_R0_helio,sc_V0_helio] = ...
%       MGALT_convertLVLHHelio(CONST,sc_R0_lvlh,sc_V0_lvlh,sc_T0_ECI_lvlh)
%
% |-----------------------------------------------------------------------
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

errorPathDisplay()
fprintf(2,'This function is not currently in use.\n');
return



% Constants
AU = CONST.AU;        	% [km/AU]
TU = CONST.TU*86400;	% [sec/TU]


% Get the coords in Helio
sc_V0 = sc_T0_ECI_lvlh'\sc_V0_lvlh;

% How do I get sc_R0???
% Need to somehow undo cross products and get norm values of stuff...

% ERROR
sc_R0 = sc_T0_ECI_lvlh'\sc_R0_lvlh;     % This method doesn't work
% ERROR

% Convert into cartesian coords from polar
% [TBD]

% R and V in DU and TU
sc_R0_helio = sc_R0 / (1/AU); 	% [AU]
sc_V0_helio = sc_V0 / (TU/AU);	% [TU]



% sc_R0_lvlh = sc_T0_ECI_lvlh(:,1);
% sc_S0_lvlh = sc_T0_ECI_lvlh(:,2);
% sc_W0_lvlh = sc_T0_ECI_lvlh(:,3);

% sc_R0_lvlh = sc_R0/sc_pos_rad;          % [DU]
% sc_W0_lvlh = cross(sc_R0,sc_V0)/norm(cross(sc_R0,sc_V0));
% sc_S0_lvlh = cross(sc_W0_lvlh,sc_R0_lvlh);
% sc_T0_ECI_lvlh = [sc_R0_lvlh, sc_S0_lvlh, sc_W0_lvlh];
% sc_V0_lvlh = sc_T0_ECI_lvlh'*sc_V0;     % [DU/TU]

% sc_pos_rad = 0;
% matrix = [0,-sc_R0_lvlh(3),sc_R0_lvlh(2);sc_R0_lvlh(3),0,-sc_R0_lvlh(1);-sc_R0_lvlh(2),sc_R0_lvlh(1),0];
% sc_R0 = sc_R0_lvlh
% *sc_pos_rad;

% % R and V in DU and TU
% sc_R0_helio = sc_R0/(1/AU);     % [DU]
% sc_V0_helio = sc_V0/(TU/AU);    % [TU]



end


