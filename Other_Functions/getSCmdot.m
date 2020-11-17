function [d_mass] = getSCmdot(CONST,OPT,T)
% FORM: [d_mass] = getSCmdot(CONST,OPT,T)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Function to get the change in spacecraft mass. Used in the 
% |     EOM ODE's
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
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
% |     -T                  (1,1)       [float]         [N]
% |         Spacecraft Thrust
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -d_mass             (1,1)       [float]         [kg/TU]
% |         The change in spacecraft mass after a timestep
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Get SC Mdot

switch OPT.thrust.mdot_method
    
    case {'constant'}       % Constant Mdot
        
        switch OPT.solver
            
            % Check if the direct method is thrusting or not
            case {'LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
                
                if T ~= 0
                    d_mass = -OPT.thrust.mdot;	% kg/s
                else
                    d_mass = 0;
                end
                
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                
                % If mdot is constant
                d_mass = -OPT.thrust.mdot;  	% kg/s
                
            otherwise
                
                errorPathDisplay()
                errorSolver()
                return
                
        end
        
    case {'equation'}
        
        % Use the already defined Solar Sail Model
        d_mass = -feval(OPT.thrust.mdot,CONST,T,OPT.thrust.Isp);  % kg/s
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect mdot_method selected.\n\n')
        return
        
end

d_mass = d_mass*(CONST.TU*86400);     % kg/TU;



end


