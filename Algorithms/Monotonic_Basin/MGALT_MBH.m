function [eval_info,selected] = ...
    MGALT_MBH(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_MBH)
% FORM: [eval_info,selected] = ...
%       MGALT_MBH(BOD,CONST,OPT,VAR,selected,num_mig,num_isl,count_MBH)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function is the Multiple Gravity-Assist Low-Thrust (MGALT) 
% |     adaptation of the Monotonic Basin Hopping (MBH) function. 
% |
% |     -This function is the main wrapper for the MBH island. Every 
% |     instance of the MBH Island object is used as an input.
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
% |     -count_MBH          (1,1)       [int]           [unitless]
% |         The current MBH island number
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



%% Pre-Allocation - Loop 1

% Solution Stagnation
stagnation_first = 0;

% Make array of 0's for the best solution and cost
optimal_soln_first = zeros(OPT.MBH(count_MBH).N1_Outer,size(VAR.bin,2));
f_best_first = 999999999999999999999*ones(OPT.MBH(count_MBH).N1_Outer,1);



%% Generate Initial Population - Loop 1

[popn_first,~] = genPopn('MBH',OPT,OPT.MBH(count_MBH),VAR);



%% Add Shared Solutions if Migration Has Occurred

if num_mig > 1
   popn_first = ISL_modReplacement(OPT.island,selected,num_mig-1,num_isl,popn_first);
end



%% First Loop

% Progress Display
clc
fprintf('~~~~~~~~~~ Migration %1.0f/%1.0f ~~~~~~~~~~\n',num_mig-1,OPT.island.Nmig)
fprintf('    Island: %1.0f/%1.0f\n',num_isl,OPT.island.Nisl);
fprintf(' Algorithm: %s\n',char(OPT.island.isl_list(num_isl,:)));
fprintf('   Loop #1: Searching Global\n');

% Pass into function
[f_first,popn_first,is_viable_first,nfeval_first] = ...
    MGALT_MBH_function(BOD,CONST,OPT,OPT.MBH(count_MBH),VAR,popn_first,num_mig,count_MBH,1);

% Generate Statistics for Current Generation
avgcost_first = mean(f_first);
mincost_first = min(f_first);
maxcost_first = max(f_first);



%% If Clusters Were Detected

% Find which rows are viable
[rows,~] = find(is_viable_first == 1);

%If no clusters were found
if isempty(rows)
    
    % Organize Best Solutions So Far
    [f_best_first,optimal_soln_first,stagnation_first] = ISL_bestSoFarCheck(f_first,popn_first,...
        f_best_first,optimal_soln_first,stagnation_first);
    
    % Selecting Solutions For Sharing
    [selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best_first,...
        optimal_soln_first,selected,num_isl,num_mig,OPT.island);
    
    % Eval Info
    eval_info.optimal_soln = sorted_optimal_soln ;
    eval_info.f_best = sorted_f_best;
    eval_info.iterations = OPT.MBH(count_MBH).N1_Outer;
    eval_info.maxcost = maxcost_first; 
    eval_info.mincost = mincost_first; 
    eval_info.avgcost = avgcost_first;
    eval_info.total_evals = nfeval_first;
    
    return
    
end    



%% Generate Initial Population - Loop 2

popn_second = MGALT_MBH_cluster(BOD,OPT,OPT.MBH(count_MBH),VAR,...
    [f_first,popn_first],rows,num_mig);



%% Pre-Allocation - Loop 2

% Solution Stagnation
stagnation_second = 0;


% Make array of 0's for the best solution and cost
optimal_soln_second = zeros(size(popn_second,1),size(popn_second,2));
f_best_second = 999999999999999999999*ones(size(popn_second,1),1);



%% Second Loop

% Progress Display
clc
fprintf('~~~~~~~~~~ Migration %1.0f/%1.0f ~~~~~~~~~~\n',num_mig-1,OPT.island.Nmig)
fprintf('    Island: %1.0f/%1.0f\n',num_isl,OPT.island.Nisl);
fprintf(' Algorithm: %s\n',char(OPT.island.isl_list(num_isl,:)));
fprintf('   Loop #2: Searching Local Basins\n');

% Pass into function
[f_second,popn_second,~,nfeval_second] = ...
    MGALT_MBH_function(BOD,CONST,OPT,OPT.MBH(count_MBH),VAR,popn_second,num_mig,count_MBH,2);
  
% Generate Statistics for Current Generation
avgcost_second = mean(f_second);
mincost_second = min(f_second);
maxcost_second = max(f_second);



%% Solutions

% Organize Best So Far Solutions
[f_best_final,optimal_soln_final,stagnation_second] = ISL_bestSoFarCheck([f_first;f_second],[popn_first;popn_second],...
    [f_best_first;f_best_second],[optimal_soln_first;optimal_soln_second],stagnation_second);

% Selecting Solutions For Sharing
[selected,sorted_optimal_soln,sorted_f_best] = ISL_modSelection(f_best_final,optimal_soln_final,selected,num_isl,num_mig,OPT.island);

% Eval Info
eval_info.optimal_soln = sorted_optimal_soln;
eval_info.f_best = sorted_f_best;
eval_info.iterations = OPT.MBH(count_MBH).N1_Outer+size(popn_second,1);
eval_info.maxcost = max(maxcost_first,maxcost_second); 
eval_info.mincost = min(mincost_first,mincost_second); 
eval_info.avgcost = mean([avgcost_first,avgcost_second]);
eval_info.total_evals = nfeval_first+nfeval_second;



end


