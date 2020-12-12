function [f_best_return,x_star_best_return] = MGALT_MBH_basin(BOD,CONST,OPT,...
    OPT_algo,VAR,f_current,x_current,num_mig,count_MBH,loop)
% FORM: [f_best_return,x_star_best_return] = MGALT_MBH_basin(BOD,CONST,OPT,...
%       OPT_algo,VAR,f_current,x_current,num_mig,count_MBH,loop)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function works alongside the "MGALT_MBH_function" to 
% |     search for basins if feasibility is detected. 
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
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         MBH option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersMBH.m"
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -f_current      	(1,1)       [float]         [unitless]
% |         The respective cost of 'x_current'
% |     -x_current       	(1,Nvar)    [float]         [unitless]
% |         The initial population which will be evaluated by the solver
% |     -num_mig            (1,1)       [int]           [unitless]
% |         The current migration number
% |     -count_MBH          (1,1)       [int]           [unitless]
% |         The current MBH island number
% |     -loop               (1,1)       [int]           [unitless]
% |         The current loop number
% |         This is if the function is solving an initial set of members 
% |         or a secondary set of members generated from MGALT_MBH_cluster
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -f_best_return    	(1,1)       [float]         [unitless]
% |         The best respective cost all x_star explored in the basin
% |     -x_star_best_return	(1,Nvar)    [float]         [unitless]
% |         The best respective member of all x_star explored in the basin
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initialize

% Number of loop
per_feas = OPT_algo.per_feas(num_mig);	% Feasibility percentage for SOI
per_rand = OPT_algo.per_rand(num_mig);	% Random percentage for generating new results
if loop == 1
    num_iter_basin = OPT_algo(count_MBH).N1_Inner;     	% Number of iterations for the inner loop
elseif loop == 2
    num_iter_basin = OPT_algo(count_MBH).N2_Outer;   	% Number of iterations for the inner loop
    per_rand = per_rand/2;
end

% ***1.0 Pre-allocate/Random point***
f_best = f_current;
x_star_best = x_current;

N_not_improve = 0;
neval = 0;



%% Search Basin

% ***4.0 If x_star is feasible***
while N_not_improve < num_iter_basin

    neval = neval+1;

    % ***4.1 Generate x_prime***
    % Generate x_prime by randomly perturbing x_current
    [x_prime] = MGALT_MBH_randomize(BOD,OPT,VAR,x_current,per_rand);

    % ***4.2 Run NLP to find x_star***
    % Run the NLP problem solver to find x_star using the initial guess x0
    [f_prime,plot_vars_prime] = feval(OPT.solver,x_prime,BOD,CONST,OPT,VAR);
    ll_Basin = fprintf('      Test: %1.0f/%1.0f\n',N_not_improve+1,num_iter_basin);

    % ***4.3 Check x_star feasibility***
    % If x_star is a feasible point AND the cost for x_star is lower than x_current
    [feasible,~,~] = MGALT_MBH_isFeasible(BOD,CONST,OPT,OPT_algo,VAR,plot_vars_prime,per_feas);
    if OPT_algo.feas_check(feasible)

        % The nested "if" loop is used to prevent the error "Operands to the 
        % || and && operators must be convertible to logical scalar values."
        if (f_prime < f_current)

            % The current stuff is now the perturbed stuff
            f_current = f_prime;
            x_current = x_prime;

            % Add onto the array of the best members and costs
            x_star_best = [x_star_best; x_prime];
            f_best = [f_best; f_prime];

            % Reset the not_improve conditions and change the percent
            % modifiers
            N_not_improve = 0;
            per_feas = per_feas/2;
            per_rand = per_rand/2;

        else
            % Add onto the not improve
            N_not_improve = N_not_improve+1;
        end

    else
        % Add onto the not improve
        N_not_improve = N_not_improve+1;
    end

    % Remove the "Test #/#"
    fprintf(repmat('\b',1,ll_Basin))

end


% ***5.0 Find the problem solution***
% Get the index of the best solution
[~,index] = min(f_best);

% Use the index to pull out the values
f_best_return = f_best(index);
x_star_best_return = x_star_best(index,:);



end


