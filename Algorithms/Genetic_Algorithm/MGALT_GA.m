function [eval_info,selected] = ...
    MGALT_GA(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_GA)
% FORM: [eval_info,selected] = ...
%       MGALT_GA(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_GA)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
% |     adaptation of the Genetic Algorithm (GA) function. 
% |
% |     -This function is the main wrapper for the GA island. Every 
% |     instance of the GA Island object is used as an input.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
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
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -selected           (Nmig,Nisl) [struct]        [unitless]
% |         Shared solutions from each island for each migration
% |     -num_mig            (1,1)       [int]           [unitless]
% |         The current migration number
% |     -num_isl            (1,1)       [int]           [unitless]
% |         The current island number
% |     -count_GA           (1,1)       [int]           [unitless]
% |         The current DE island number
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -eval_info          (1,1)       [struct]        [unitless]
% |         Best solutions, best costs, number of iterations, etc...
% |     -selected           (Nmig,Nisl) [struct]        [unitless]
% |         Shared solutions from each island for each migration
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Pre-Allocation

% Cost Parameters
f = zeros(OPT.GA(count_GA).Npop,1);
avgcost = zeros(OPT.GA(count_GA).Ngen,1);
mincost = zeros(OPT.GA(count_GA).Ngen,1);
maxcost = zeros(OPT.GA(count_GA).Ngen,1);

% Solution Stagnation
stagnation = 0; 

% Number of evals
nfeval = 0;

% Make arrays of 0's for the best solution and cost
optimal_soln = zeros(OPT.GA(count_GA).Npop,size(VAR.bin,2));
f_best = 999999999999999999999*ones(OPT.GA(count_GA).Npop,1);



%% Generate Initial Population

[popn,~] = genPopn('GA',OPT,OPT.GA(count_GA),VAR);



%% Add Shared Solutions if Migration Has Occurred

if num_mig > 1
   popn = ISL_modReplacement(OPT.island,selected,num_mig-1,num_isl,popn);
end



%% Iterate Generations

% Progress Display
clc
fprintf('~~~~~~~~~~~ Migration %1.0f/%1.0f ~~~~~~~~~~~\n',num_mig-1,OPT.island.Nmig)
fprintf('    Island: %1.0f/%1.0f\n',num_isl,OPT.island.Nisl);
fprintf(' Algorithm: %s\n',char(OPT.island.isl_list(num_isl,:)));

for num_gen = 1:OPT.GA(count_GA).Ngen
    
    % Calculate Fitness Value for Each Member
    if OPT.parallel
        
        % Preallocate
        data(1:OPT.GA(count_GA).Npop) = parallel.FevalFuture;
        fh = str2func(string(OPT.solver));      % https://www.mathworks.com/help/matlab/ref/str2func.html
        
        % Display Info
        ll_par = fprintf('\n    Parallel Processing   \n');
        ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n');
        ll_gen = fprintf('Generation: %1.0f/%1.0f\n',num_gen,OPT.GA(count_GA).Ngen);
        
        % Cost
        for member = 1:OPT.GA(count_GA).Npop
            data(member) = parfeval(OPT.ppool,fh,1,popn(member,:),BOD,CONST,OPT,VAR);
        end
        
        % Collect results as they become available
        for member = 1:OPT.GA(count_GA).Npop
            [index,cost] = fetchNext(data);
            f(index) = cost;
        end
        fprintf(repmat('\b',1,ll_gen))
        fprintf(repmat('\b',1,ll_RAM))
        fprintf(repmat('\b',1,ll_par))
        
    else
    
        % Display Info
        ll_gen = fprintf('Generation: %1.0f/%1.0f\n',num_gen,OPT.GA(count_GA).Ngen);
            
        for member = 1:OPT.GA(count_GA).Npop
            % Cost
            [f(member),~] = feval(OPT.solver,popn(member,:),BOD,CONST,OPT,VAR);
        end
        fprintf(repmat('\b',1,ll_gen))
    
    end
    nfeval = nfeval + OPT.GA(count_GA).Npop;
    
    % Generate Statistics for Current Generation
    avgcost(num_gen) = mean(f); 
    mincost(num_gen) = min(f); 
    maxcost(num_gen) = max(f);
  
    % Organize Best So Far Solutions
    [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation);
  
    % Select Mating Pool; Mate to Create Next Generation
    popn = MGALT_GA_nextGeneration(BOD,OPT,OPT.GA(count_GA),VAR,popn,f);
    
end



%% Selecting Solutions For Sharing

[selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best,optimal_soln,selected,num_isl,num_mig,OPT.island);



%% Eval Info

eval_info.optimal_soln = sorted_optimal_soln;
eval_info.f_best = sorted_f_best;
eval_info.iterations = OPT.GA(count_GA).Ngen;
eval_info.maxcost = maxcost; 
eval_info.mincost = mincost; 
eval_info.avgcost = avgcost;
eval_info.total_evals = nfeval;



end


