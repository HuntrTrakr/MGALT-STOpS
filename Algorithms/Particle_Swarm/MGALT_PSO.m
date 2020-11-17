function [eval_info,selected] = ...
    MGALT_PSO(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_PSO)
% FORM: [eval_info,selected] = ...
%       MGALT_PSO(solver,inputs,opt_PSO,selected,mig,isl)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
% |     adaptation of the Particle Swarm Optimization (PSO) function. 
% |
% |     -This function is the main wrapper for the PSO island. Every 
% |     instance of the PSO Island object is used as an input.
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
% |     -count_PSO          (1,1)       [int]           [unitless]
% |         The current PSO island number
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
f = zeros(OPT.PSO(count_PSO).Npop,1);
avgcost = zeros(OPT.PSO(count_PSO).Npop,1);
mincost = zeros(OPT.PSO(count_PSO).Npop,1);
maxcost = zeros(OPT.PSO(count_PSO).Npop,1);

% Solution Stagnation
stagnation = 0; 

% Number of evals
nfeval = 0;

% Make arrays of 0's for the best solution and cost
optimal_soln = zeros(OPT.PSO(count_PSO).Npop,size(VAR.bin,2));
f_best = 999999999999999999999*ones(OPT.PSO(count_PSO).Npop,1);



%% Generate Initial Population

[popn,bee] = genPopn('PSO',OPT,OPT.PSO(count_PSO),VAR);



%% Add Shared Solutions if Migration Has Occurred

if num_mig > 1
   popn = ISL_modReplacement(OPT.island,selected,num_mig-1,num_isl,popn);
   for p = 1:OPT.PSO(count_PSO).Npop
      bee(p).pos = popn(p,:);
   end
end



%% Iterate Generations

% Progress Display
clc
fprintf('~~~~~~~~~~~ Migration %1.0f/%1.0f ~~~~~~~~~~~\n',num_mig-1,OPT.island.Nmig)
fprintf('    Island: %1.0f/%1.0f\n',num_isl,OPT.island.Nisl);
fprintf(' Algorithm: %s\n',char(OPT.island.isl_list(num_isl,:)));
    
for t = 1:OPT.PSO(count_PSO).tspan
    
    % Calculate Fitness Value for Each Member
    if OPT.parallel
    
        % Preallocate
        data(1:OPT.PSO(count_PSO).Npop) = parallel.FevalFuture;
        fh = str2func(string(OPT.solver));      % https://www.mathworks.com/help/matlab/ref/str2func.html
        
        % Display Info
        ll_par = fprintf('\n    Parallel Processing   \n');
        ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n');
        ll_time = fprintf(' Time Step: %1.0f/%1.0f\n',t,OPT.PSO(count_PSO).tspan); 
        
        % Cost
        for member = 1:OPT.PSO(count_PSO).Npop
            data(member) = parfeval(OPT.ppool,fh,1,popn(member,:),BOD,CONST,OPT,VAR);
        end
        
        % Collect results as they become available
        for member = 1:OPT.PSO(count_PSO).Npop
            [index,cost] = fetchNext(data);
            bee(index).f = cost;
            f(index) = cost;
        end
        fprintf(repmat('\b',1,ll_time))
        fprintf(repmat('\b',1,ll_RAM))
        fprintf(repmat('\b',1,ll_par))

    else
    
        % Display Info
        ll_time = fprintf(' Time Step: %1.0f/%1.0f\n',t,OPT.PSO(count_PSO).tspan);
        
        for member = 1:OPT.PSO(count_PSO).Npop
           [bee(member).f, ~] = feval(OPT.solver,bee(member).pos,BOD,CONST,OPT,VAR); 
           f(member) = bee(member).f;
        end
        fprintf(repmat('\b',1,ll_time))
    
    end
    
    % See If Fitness Is Best Yet
    for member = 1:OPT.PSO(count_PSO).Npop
        if f(member) < bee(member).f_p
        	bee(member).p   = bee(member).pos;
            bee(member).f_p = f(member);
        end
    end
    
    % Generate Statistics for Current Time
    avgcost(t) = mean(f);
    mincost(t) = min(f);
    maxcost(t) = max(f);
    
    % Organize Best So Far Solutions
    [f_best,optimal_soln,stagnation] = ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation);
    
    % Population Members Move, Communicate, & Adjust Velocities
    bee = MGALT_PSO_nextGeneration(BOD,OPT,OPT.PSO(count_PSO),VAR,bee);
    
    for p = 1:OPT.PSO(count_PSO).Npop
        popn(p,:) = bee(p).pos;
    end
    
    nfeval = nfeval + OPT.PSO(count_PSO).Npop;
    
end



%% Select Solutions For Sharing

[selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best,optimal_soln,selected,num_isl,num_mig,OPT.island);



%% Eval Info

eval_info.optimal_soln = sorted_optimal_soln;
eval_info.f_best = sorted_f_best;
eval_info.iterations = OPT.PSO(count_PSO).tspan;
eval_info.maxcost = maxcost; 
eval_info.mincost = mincost; 
eval_info.avgcost = avgcost;
eval_info.total_evals = nfeval;



end


