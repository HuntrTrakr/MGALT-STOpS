function [x_prime] = MGALT_MBH_randomize(BOD,OPT,VAR,x_current,per_rand)
% FORM: [x_prime] = MGALT_MBH_randomize(BOD,OPT,VAR,x_current,per_rand)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function creates random perturbations within x_current and
% |     returns it as x_prime
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -x_current          (1,Nvar)    [float]         [unitless]
% |         The current population which will be perturbed
% |     -per_rand           (1,1)       [float]     	[unitless]
% |         Random percentage number for generating new results
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -x_prime            (1,Nvar)    [float]         [unitless]
% |         The final perturbed population
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Perform Randomization

% Perallocate
x_prime = zeros(1,size(x_current,2));

% Randomization
switch OPT.solver
        
    case {'MGALT_DIR_FSM_2D','MGALT_DIR_FBSM_2D'}
        
        % Apply perturbations on every member
        for i1 = 1:length(x_prime)
            random_gen = x_current(i1)*per_rand;
            
            if VAR.bin(i1)   % Binary variables
                x_prime(i1) = randomNum(x_current(i1)-random_gen,...
                x_current(i1)+random_gen, 'int');
            else
                x_prime(i1) = randomNum(x_current(i1)-random_gen,...
                x_current(i1)+random_gen, 'dec');
            end
            
        end
        
        % Ensure the first thrust = 1 or 0
        switch OPT.thrust.thrust_method
            
            case {'constant_thrust','equation_thrust'}
                
                x_prime(2) = 1;     % Ensure that always thrusting when leaving departure planet
                
            case {'variable_thrust'}
                
                % Acts as a pass
                
            otherwise
                
                errorPathDisplay();
                fprintf(2,'Incorrect thrust method selected.\n');
                return
                
        end
   
    case {'MGALT_IN_FSM_2D','MGALT_IN_FBSM_2D'}

        % Apply perturbations on every member
        for i1 = 1:length(x_prime)
            random_gen = x_current(i1)*per_rand;
            x_prime(i1) = randomNum(x_current(i1)-random_gen,...
                x_current(i1)+random_gen, 'dec');
        end
        
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end

% Ensure the bounds of the search space haven't been exceeded
for i5 = 1:length(x_prime)
    
    if x_prime(i5) < VAR.low(i5)
        x_prime(i5) = VAR.low(i5);
    end
    
    if x_prime(i5) > VAR.high(i5)
        x_prime(i5) = VAR.high(i5);
    end
    
end


% Ensure that ToF is not broken
x_prime = MGALT_fixToF(BOD,OPT,VAR,x_prime);



end


