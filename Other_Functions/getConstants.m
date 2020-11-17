function [CONST] = getConstants(BOD)
% FORM: [CONST] = getConstants(BOD)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Get the constants for this problem
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -CONST              (1,1)       [struct]        [unitless]
% |         A struct containing constants used in the calcs. Contains
% |         values for AU, TU, Sun (rad/mu/rp) and (rad/mu/rp/SOI/per) 
% |         for any bodies used in the optimization scheme. This is a 
% |         dynamic struct and will adapt to contain only the necesary 
% |         information
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Get Constants

bodies = BOD.bodies;

% Get the canonical units
CONST.AU = constants('AU');
CONST.TU = constants('TU');

% Get sun units
CONST.Sun_rad = constants('radius Sun');
CONST.Sun_mu = constants('mu Sun');
CONST.Sun_rp = constants('min rp Sun');


% Get all related info to bodies
for i1 = 1:length(bodies)
    
    switch bodies{i1}
        
        case {'Mercury'}
            
            CONST.Mercury_rad = constants('radius Mercury');
            CONST.Mercury_mu = constants('mu Mercury');
            CONST.Mercury_rp = constants('min rp Mercury');
            CONST.Mercury_SOI = constants('SOI Mercury');
            CONST.Mercury_per = constants('period Mercury');
            
        case {'Venus'}
            
            CONST.Venus_rad = constants('radius Venus');
            CONST.Venus_mu = constants('mu Venus');
            CONST.Venus_rp = constants('min rp Venus');
            CONST.Venus_SOI = constants('SOI Venus');
            CONST.Venus_per = constants('period Venus');
            
        case {'Earth'}
            
            CONST.Earth_rad = constants('radius Earth');
            CONST.Earth_mu = constants('mu Earth');
            CONST.Earth_rp = constants('min rp Earth');
            CONST.Earth_SOI = constants('SOI Earth');
            CONST.Earth_per = constants('period Earth');
            
        case {'Mars'}
            
            CONST.Mars_rad = constants('radius Mars');
            CONST.Mars_mu = constants('mu Mars');
            CONST.Mars_rp = constants('min rp Mars');
            CONST.Mars_SOI = constants('SOI Mars');
            CONST.Mars_per = constants('period Mars');
            
        case {'Jupiter'}
            
            CONST.Jupiter_rad = constants('radius Jupiter');
            CONST.Jupiter_mu = constants('mu Jupiter');
            CONST.Jupiter_rp = constants('min rp Jupiter');
            CONST.Jupiter_SOI = constants('SOI Jupiter');
            CONST.Jupiter_per = constants('period Jupiter');
            
        case {'Saturn'}
            
            CONST.Saturn_rad = constants('radius Saturn');
            CONST.Saturn_mu = constants('mu Saturn');
            CONST.Saturn_rp = constants('min rp Saturn');
            CONST.Saturn_SOI = constants('SOI Saturn');
            CONST.Saturn_per = constants('period Saturn');
            
        case {'Uranus'}
            
            CONST.Uranus_rad = constants('radius Uranus');
            CONST.Uranus_mu = constants('mu Uranus');
            CONST.Uranus_rp = constants('min rp Uranus');
            CONST.Uranus_SOI = constants('SOI Uranus');
            CONST.Uranus_per = constants('period Uranus');
            
        case {'Neptune'}
            
            CONST.Neptune_rad = constants('radius Neptune');
            CONST.Neptune_mu = constants('mu Neptune');
            CONST.Neptune_rp = constants('min rp Neptune');
            CONST.Neptune_SOI = constants('SOI Neptune');
            CONST.Neptune_per = constants('period Neptune');
            
        case {'Pluto'}
            
            CONST.Pluto_rad = constants('radius Pluto');
            CONST.Pluto_mu = constants('mu Pluto');
            CONST.Pluto_rp = constants('min rp Pluto');
            CONST.Pluto_SOI = constants('SOI Pluto');
            CONST.Pluto_per = constants('period Pluto');
            
        otherwise
            
            errorPathDisplay();
            fprintf(2,'Incorrect body selection.\n')
            return
            
    end
    
end



end


