function [T] = getSCThrust(CONST,OPT,pos_rad,mem_thrust)
% FORM: [T] = getSCThrust(CONST,OPT,pos_rad,mem_thrust)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function to get the current thrust value of the spacecraft. Used
% |     in the EOM ODE's
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
% |     -pos_rad            (1,1)       [float]         [AU]
% |         The radial position of the S/C with respect to the sun
% |     -mem_thrust         (1,n)       [boolean/float] [N]
% |         Optional, for the direct method to let STOpS know if there is
% |         thrust during that time segment (0 or 1), or if using variable 
% |         thrust, the amount of thrust
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -T                  (1,1)       [float]         [N]
% |         Spacecraft Thrust
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Get SC Thrust

switch OPT.thrust.thrust_method
    
    case {'constant'}       % Constant Thrust
        
        switch OPT.solver
            
            % Check if the direct method is thrusting or not
            case {'LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
                
                if mem_thrust
                    T = OPT.thrust.thrust;	% N
                else
                    T = 0;                  % N
                end
                
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                
                T = OPT.thrust.thrust;      % N

            otherwise
                
                errorPathDisplay()
                errorSolver()
                return
                
        end
        
    case {'variable'}
        
        switch OPT.solver
            
            case {'LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
                
                T = mem_thrust;     % N
                
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                
                errorPathDisplay()
                fprintf(2,'The INDIRECT Method does not support varaible thrust.\n')
                return
                
            otherwise
                
                errorPathDisplay()
                errorSolver()
                return
                
        end
        
    case {'equation'}
        
        switch OPT.solver
            
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                
                T = feval(OPT.thrust.thrust,OPT,pos_rad);       % N
                
            case {'LT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
                
                if mem_thrust
                    T = feval(OPT.thrust.thrust,OPT,pos_rad);	% N
                else
                    T = 0;                                      % N
                end
                
            otherwise
                
                errorPathDisplay()
                errorSolver()
                return
                
        end
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect thrust profile selected.\n\n')
        return
        
end

T = T/1000;         % (kg*km/s^2)
T = T*(CONST.TU*86400)^2/CONST.AU;      % (kg*DU/TU^2)



end


