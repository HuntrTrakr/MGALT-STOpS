function [popn_final] = ...
    MGALT_MBH_cluster(BOD,OPT,OPT_algo,VAR,popn_initial,rows,num_mig)
% FORM: [popn_final] = ...
%       MGALT_MBH_cluster(BOD,OPT,OPT_algo,VAR,popn_initial,rows,num_mig)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function works alongside the "MGALT_MBH" function to 
% |     generate new population members from possible clusters found 
% |     during the first iteration of MBH.
% |
% |     -If the first loop determines that there are clusters from the 
% |     "MGALT_MBH_isFeasible" function, those members of the population 
% |     are used and mutated to form a set of new members for a new 
% |     population subset
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
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         MBH option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersMBH.m"
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -popn_initial	(N1_outer,Nvar+1)   [float] 	[unitless]
% |         The initial population which was run during the first loop, 
% |         bundled with the respective cost of the population members
% |     -rows               (1,n)       [int]           [unitless]
% |         The row numbers where feasible solutions are
% |     -num_mig            (1,1)       [int]           [unitless]
% |         The current migration number
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -popn_final     (N2_outer*clst,Nvar) [float]    [unitless]
% |         The final population which was generated from the best results 
% |         of the popn_initial. This contains a lot more members beause 
% |         each member of the initial was mutated N2_outer times
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Unpack Data

num_iter = OPT_algo.N2_Outer;
max_clst = OPT_algo.maxclst;
per_rand = OPT_algo.per_rand(num_mig);



%% Function

% ***1.0 Data in the cluster***

% Data from the rows which were viable
data_indexed = popn_initial(rows,:);

% Sort by lowest cost
[~,idx] = sort(data_indexed(:,1));
sorted = data_indexed(idx,:);
bounds_second = sorted(:,2:end);

% Adjust for the number of clusters
if size(bounds_second,1) <= max_clst
    max_clst = size(bounds_second,1);
end


% ***2.0 Generate Population***
popn_final = zeros(num_iter*max_clst,size(bounds_second,2));

% Split the initial population into different chunks 
for i1 = 1:max_clst
    for i2 = 1:num_iter
        popn_final(i2+((i1-1)*num_iter),:) = MGALT_MBH_randomize(BOD,OPT,...
            VAR,bounds_second(i1,:),per_rand/2);
    end
end
      


end


