function [opt_cost,opt_weight,opt_ode] = optionsCost(OPT,selection)
% FORM: [opt_cost,opt_weight,opt_ode] = optionsCost(OPT,selection)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -These paramaters dictate how the cost will be calculated. 
% |     The goal of the cost function is to have the final condtions for 
% |     a member be as close to the desired final condtions as possible 
% |     and to minimize the time of the trajectory.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -selection          (1,n)       [string]        [unitless]
% |         A used defined string which is used in a switch/case format
% |         to return predefined options
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -opt_cost           (1,n)       [assorted]      [unitless]
% |         The return object with cost parameters
% |     -opt_weight         (1,n)       [assorted]      [unitless]
% |         The return object with weighting parameters
% |     -opt_ode            (1,n)       [assorted]      [unitless]
% |         The return object with ode parameters
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     The parameters are defined as follows:
% |
% |
% |     -tolR:              (1,1)       [float]         [percent]
% |         Radius convergence factor from 0-1. A value of 0.01 is 
% |         recommended because that means the cost is driving to have 
% |         the radius converge within 1% of the desired value.
% |     -tolTheta:          (1,1)       [float]         [percent]
% |         Theta convergence factor from 0-1. A value of 0.1 is 
% |         recommended because that means the cost is driving to have 
% |         the angular displacement converge within 10% of the desired 
% |         value.
% |     -tolU:              (1,1)       [float]         [percent]
% |         Radial velocity convergence factor from 0-1. A value of 0.01 
% |         is recommended because that means the cost is driving to have 
% |         the radial velocity converge within 1% of the desired value.
% |     -tolV:              (1,1)       [float]         [percent]
% |         Tangential velocity convergence factor from 0-1. A value of 
% |         0.01 is recommended because that means the cost is driving to 
% |         have the tangential velocity converge within 1% of the desired 
% |         value.
% |     -R:                 (1,1)       [int]           [binary]
% |         Radial cost switch. Choose 1 to include radius convergence in 
% |         the cost calculation.
% |     -Theta:             (1,1)       [imt]           [binary]
% |         Theta cost switch. Choose 1 to include angular displacement 
% |         convergence in the cost calculation. Choosing 1 result in the 
% |         trajectory terminating at the arrival planet while choosing 0 
% |         will result in the trajectory terminating anywhere on the 
% |         orbit of the the arrival planet. Choose 0 if an orbit 
% |         transfer is desired.
% |     -U:                 (1,1)       [int]           [binary]
% |         Radial velocity cost switch. Choose 1 to include radial 
% |         velocity convergence in the cost calculation.
% |     -V:                 (1,1)       [int]           [binary]
% |         Tangential velocity cost switch. Choose 1 to include 
% |         tangential velocity convergence in the cost calculation.
% |     -tt:                (1,1)       [int]           [binary]
% |         End time cost switch. Choose 1 to include time convergence 
% |         in the cost calculation.
% |     -W1:                (1,1)       [float]         [unitless]
% |         Weight for end time convergence. Choosing a bigger number 
% |         allows the end time to vary more. Choosing a smaller value 
% |         puts more pressure on minimizing the time, but risks 
% |         sacrificing the convergence on the desired position and 
% |         velocity. The default value is 3.5.
% |
% |-----------------------------------------------------------------------



%% Selection

opt_cost = struct;
opt_weight = struct;
opt_ode = struct;

switch OPT.solver
    
    case {'LT_DIR_FSM_2D','LT_IN_FSM_2D'}
        
        switch selection
            
            case {'default'}
                
                % Cost
                opt_cost.tolR = 0.01;           % Radius convergence convergence tolerance
                opt_cost.tolTheta = 0.1;        % Angular displacement convergence tolerance
                opt_cost.tolU = 0.01;           % Radial velocity convergence tolerance
                opt_cost.tolV = 0.01;           % Tangential velocity convergence tolerance
                opt_cost.R = 1;                 % Radial cost switch, on or off
                opt_cost.Theta = 1;             % Angular displacement cost switch, on or off
                opt_cost.U = 1;                 % Radial velocity cost switch, on or off
                opt_cost.V = 1;                 % Tangential velocity cost switch, on or off
                opt_cost.tt = 1;                % End time cost switch, on or off

                % Weighting
                opt_weight.W_tof_conv = 3.5;	% Weight factor for end time convergence
                
                % ODE
                opt_ode = odeset('AbsTol',1e-6,'RelTol',1e-6);
                
            case {'future_cases'}
                
                fprintf('Add in future choices here\n')
                
            otherwise
                
                errorPathDisplay();
                errorSolver();
                return
                
        end
        
    case {'MGALT_DIR_FBSM_2D','MGALT_IN_FBSM_2D'}
        
        switch selection
            
            case {'default'}
                
                % Cost
                opt_cost.tolR = 0.001;          % Radius convergence convergence tolerance
                opt_cost.tolTheta = 0.05;       % Angular displacement convergence tolerance
                opt_cost.tolU = 0.001;          % Radial velocity convergence tolerance
                opt_cost.tolV = 0.001;          % Tangential velocity convergence tolerance
                opt_cost.R = 1;                 % Radial cost switch, on or off
                opt_cost.Theta = 1;             % Angular displacement cost switch, on or off
                opt_cost.U = 1;                 % Radial velocity cost switch, on or off
                opt_cost.V = 1;                 % Tangential velocity cost switch, on or off
                opt_cost.tt = 1;                % End time cost switch, on or off
                opt_cost.dV = 1;                % dV for v_inf1 and v_inf2, on or off

                % Weighting
                opt_weight.W_tof_conv = 3.5;   	% Weight factor for end time convergence
                opt_weight.control_v = 0.005;	% Weight factor for the tolerance for the dV step
                opt_weight.MBH_feas = 1e-5;    	% Weight factor for determining if the solution falls into feasibility
                    % ^ Really should be feasibility for the whole thing
                    % for whoever works on this after me
                    
                % ODE
                opt_ode = odeset('AbsTol',1e-6,'RelTol',1e-6);
                
            case {'future_cases'}
                
                fprintf('Add in future choices here\n')
                
            otherwise
                
                errorPathDisplay();
                errorSolver();
                return
                
        end
        
end



end


