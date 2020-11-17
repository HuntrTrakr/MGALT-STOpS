function [f_best,optimal_soln,stagnation] = ...
    ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
% FORM: [f_best,optimal_soln,stagnation] = ...
%       ISL_bestSoFarCheck(f,popn,f_best,optimal_soln,stagnation)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Best So Far Check funciton. Sorts the solutions into the best
% |     order
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -f                  (Npop,1)  	[float]      	[unitless]
% |         The cost associated with every variable string in the 
% |         population
% |     -popn               (Npop,Nvar) [float]         [unitless]
% |         The current population of Npop variable strings
% |     -f_best             (Npop,1)	[float]         [unitless]
% |         The cost associated with every variable string in the 
% |         population in best to worst order
% |     -optimal_soln    	(Npop,Nvar) [float]         [unitless]
% |         New array of optimal variables strings associated with costs 
% |         in 'f_best'
% |     -stagnation       	(1,1)       [int]           [unitless]
% |         Measure of if the algorithm is stagnating. Counts how many 
% |         generations in a row with no better solution
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -f_best             (Npop,1)  	[float]      	[unitless]
% |         New array of best costs in order
% |     -optimal_soln    	(Npop,Nvar) [float]         [unitless]
% |         New array of optimal variables strings associated with costs 
% |         in 'f_best'
% |     -stagnation       	(1,1)       [int]           [unitless]
% |         Measure of if the algorithm is stagnating. Counts how many 
% |         generations in a row with no better solution
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

Npop = length(f);
found_better = 0;

for j = 1:Npop
    
   if isempty(f)
       break;
   end
   
   while any(min(f) == f_best)
      [~,cut] = min(f);
      f(cut) = []; 
      popn(cut,:) = [];
      
      if isempty(f)
          break
      end
      
   end
   
   if isempty(f)
       break
   end
   
   if any(min(f) < f_best) % If it's better than any of them
      found_better = 1;
      [~,replace_loc] = max(f_best);
      [f_best(replace_loc),ind] = min(f);
      optimal_soln(replace_loc,:) = popn(ind,:);
      f(ind) = []; 
      popn(ind,:) = [];
   end
   
end

if found_better
    stagnation = 0;
else
    stagnation = stagnation + 1; % If a better tour was not found
end



end


