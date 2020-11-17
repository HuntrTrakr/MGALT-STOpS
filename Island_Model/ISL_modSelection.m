function [selected,sorted_optimal_soln,sorted_f_best] = ...
    ISL_modSelection(f_best,optimal_soln,selected,num_isl,nmu_mig,OPT_island)
% FORM: [selected,sorted_optimal_soln,sorted_f_best] = ...
%       ISL_modSelection(f_best,optimal_soln,selected,isl,mig,opt_Isl)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Section function for solution sharing between islands. This 
% |     function selects the desired solutions from an single islands 
% |     population after that island runs
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -f_best             (Npop,1)	[float]         [unitless]
% |         The cost associated with every variable string in the 
% |         population in best to worst order
% |     -optimal_soln       (Npop,Nvar) [float]         [unitless]
% |         The current population of Npop variable strings that 
% |         coorespond with 'f_best'
% |     -selected        	(Nmig,Nisl)	[struct]        [unitless]
% |         Shared solutions from each island for each migration
% |     -num_isl            (1,1)       [int]           [unitless]
% |         The current island number
% |     -num_mig            (1,1)       [int]           [unitless]
% |         The current migration number
% |     -OPT_island     	(1,1)       [struct]        [unitless]
% |         Island Model option parameters
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -selected           (Nmig,Nisl) [struct]        [unitless]
% |         Shared solutions from each island for each migration
% |     -sorted_optimal_soln	(Npop,Nvar) [float]  	[unitless]
% |         The current population of Npop variable strings that 
% |         coorespond with 'sorted_f_best'
% |     -sorted_f_best  	(Npop,1)	[float]         [unitless]
% |         The cost associated with every variable string in the 
% |         population in best to worst order
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

[Npop,~] = size(optimal_soln);

% Make sure N is less than or equal to population size
if ~strcmpi(OPT_island.sel_pol(num_isl,:),'threshold') 
    if OPT_island.sel_opt(num_isl) > Npop
        OPT_island.sel_opt(num_isl) = Npop;
    end
end

% Sort Best to Worst
[sorted_f_best,I] = sort(f_best);
sorted_optimal_soln = optimal_soln(I,:);

% Preallocate for Selected Solutions
selected(nmu_mig,num_isl).f = [];
selected(nmu_mig,num_isl).soln = [];



%% Selection

switch char(OPT_island.sel_pol(num_isl,:))
    
    case {'random'}                         % Random
        
        Nsel = OPT_island.sel_opt(num_isl);
        rand_solns = optimal_soln; 
        rand_f = f_best;
        
        for j = 1:Nsel
            choice = randomNum(1,Npop-j+1,'int');
            selected(nmu_mig,num_isl).soln(j,:) = rand_solns(choice,:);
            selected(nmu_mig,num_isl).f(j) = rand_f(choice);
            rand_solns(choice,:) = [];
            rand_f(choice) = [];
            if isempty(rand_solns)
                break
            end
        end
        
    case {'natural_selection'}               % Natural Selection
        
        Nsel = OPT_island.sel_opt(num_isl);
        for j = 1:Nsel
            selected(nmu_mig,num_isl).soln(j,:) = sorted_optimal_soln(j,:); 
            selected(nmu_mig,num_isl).f(j) = sorted_f_best(j);
        end
        
    case {'threshold'}                       % Threshold
        
        thresh = OPT_island.sel_opt(num_isl);
        j = 1;
        while sorted_f_best(j) <= thresh && j <= Npop 
            selected(nmu_mig,num_isl).soln(j,:) = sorted_optimal_soln(j,:);
            selected(nmu_mig,num_isl).f(j) = sorted_f_best(j);
            j = j + 1;
            if j > length(sorted_f_best)
                break
            end
        end
        
    case {'rank_weighted','cost_weighted'}   % Rank or Cost Weighted
        
        if strcmpi(OPT_island.sel_pol(num_isl,:),'rank_weighted')
            den = sum(1:length(sorted_f_best));
            for i = 1:length(sorted_f_best)
                p(i) = (length(sorted_f_best)-i+1)/den;
            end
        elseif strcmpi(OPT_island.sel_pol(num_isl,:),'cost_weighted')
            den = sum(sorted_f_best);
            for i = 1:length(sorted_f_best)
                p(i) = sorted_f_best(i)/den;
            end
        end
        
        Nsel = OPT_island.sel_opt(num_isl);
        
        for k = 1:Nsel
            spin = randomNum( 0, sum(p), 'dec' );
            for j = 1:length(p)
                if spin <= sum(p(1:j))
                    selected(nmu_mig,num_isl).soln(k,:) = sorted_optimal_soln(j,:);
                    selected(nmu_mig,num_isl).f(k) = sorted_f_best(j);
                    break
                end
            end
        end
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect Island Model Selection method selected.\n')
        return
        
end



end


