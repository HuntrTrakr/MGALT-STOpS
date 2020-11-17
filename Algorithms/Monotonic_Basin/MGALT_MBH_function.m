function [f_best_return,x_star_best_return,is_viable,neval] = ...
    MGALT_MBH_function(BOD,CONST,OPT,OPT_algo,VAR,x0,num_mig,count_MBH,loop)
% FORM: [f_best_return,x_star_best_return,is_viable,neval] = ...
%       MGALT_MBH_function(BOD,CONST,OPT,OPT_algo,VAR,x0,num_mig,count_MBH,loop)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function is the actual implimentation of the MBH algorithm. 
% |     The algorithm had to be broken into a function because it would 
% |     not be feasible to fit into "MGALT_MBH" due to organization
% |
% |     -A user who would like solutions where ALL of the transfer 
% |     trajectories are feasible needs to chenge the comment under 
% |     section 3.0 and 4.3 from "any(feasible)" to "all(feasible)". The 
% |     rational for having any was explained in the accompanying thesis, 
% |     as the randomness could potentially alows MBH to find a solution 
% |     where both transfers are considered feasible. As on now, it will 
% |     accept a solution where any transfer is feasible
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
% |     -x0                 (1,Nvar)    [float]         [unitless]
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
% |     -f_best             (1,1)       [float]         [unitless]
% |         The respective cost of each member in 'x0'
% |     -member_final       (1,Nvar)    [float]         [unitless]
% |         The new perturbed member consisting of the original member and
% |         the perturbations added to it
% |     -is_viable          (1,1)       [boolean]   	[unitless]
% |         The new perturbed member consisting of the original member and
% |         the perturbations added to it
% |     -nfeval           	(1,1)       [int]           [unitless]
% |         The current number of iterations for this MBH object
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initialize

% Number of iterations
num_iter_global = size(x0,1);     % Number of iterations for the outer loop

% Number of loop
per_feas = OPT_algo.per_feas(num_mig);	% Feasibility percentage for SOI
per_rand = OPT_algo.per_rand(num_mig);	% Random percentage for generating new results
if loop == 1
    num_iter_basin = OPT_algo(count_MBH).N1_Inner;     	% Number of iterations for the inner loop
    disp_info = ' ---Potential Basin Found---';
elseif loop == 2
    num_iter_basin = OPT_algo(count_MBH).N2_Outer;   	% Number of iterations for the inner loop
    per_rand = per_rand/2;
    disp_info = ' ---Exploring Minima---';
end

% Preallocate Cost Parameters
f = zeros(num_iter_global,1);
f_best_return = zeros(num_iter_global,1);
x_star_best_return = zeros(num_iter_global,size(x0,2));
is_viable = zeros(num_iter_global,1);

% Perallocate the plot_vars struct. Need to run 1 execution to get field
% names
[~,vars] = feval(OPT.solver,x0(1,:),BOD,CONST,OPT,VAR);
name = fieldnames(vars);
for i0 = 1:length(name)
    plot_vars(num_iter_global).(name{i0}) = [];
end



%% feval
    
if OPT.parallel
    
    % Function Handle
    fh = str2func(string(OPT.solver));      % https://www.mathworks.com/help/matlab/ref/str2func.html
    
    % Display info
    ll_par = fprintf('\n    Parallel Processing   \n');
    ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n');

    % Chunks to break parallel into to prevent a lot of RAM usage at once
    chunk_nums = 1000;      % Set by Malloy, uses ~2.0GB of RAM per pass
    chunk_whole = floor(num_iter_global/chunk_nums);
    chunk_rems = mod(num_iter_global,chunk_nums);
    
    % Run through all searches minus mod
    for chunk_iter = 1:chunk_whole
        
        % Preallocate temp vars
        f_temp = zeros(chunk_nums,1);
        for i0 = 1:length(name)
            plot_vars_temp(chunk_nums).(name{i0}) = [];
        end
        data_temp(1:chunk_nums) = parallel.FevalFuture;
        
        % Display info
        ll_loop1 = fprintf('    Search: %1.0f-%1.0f / %1.0f\n',((chunk_iter-1)*chunk_nums+1),(chunk_iter*chunk_nums),num_iter_global);
        ll_loop2 = fprintf('    Search: %1.0f\n',(chunk_iter-1)*chunk_nums+1);
        
        % Perform parfeval on chunk n
        for chunk_num = 1:chunk_nums
            
            % Search Display
            fprintf(repmat('\b',1,ll_loop2))
            ll_loop2 = fprintf('    Search: %1.0f\n',((chunk_iter-1)*chunk_nums)+chunk_num); 
            
            data_temp(chunk_num) = parfeval(OPT.ppool,fh,2,x0(((chunk_iter-1)*chunk_nums)+chunk_num,:),BOD,CONST,OPT,VAR);
        end
        fprintf(repmat('\b',1,ll_loop2))
        
        % Collect results for chunk n
        ll_collect = fprintf('\n ---Collecting results from parpool---\n');
        for chunk_num = 1:chunk_nums
            [index,cost,plt_vars] = fetchNext(data_temp);
            f_temp(index) = cost;
            plot_vars_temp(index) = plt_vars;
        end
        
        % Save the collected results into the total struct
        f( ((chunk_iter-1)*chunk_nums+1):(chunk_iter*chunk_nums) ) = f_temp;
        plot_vars( ((chunk_iter-1)*chunk_nums+1):(chunk_iter*chunk_nums) ) = plot_vars_temp;
        
        % Clear display
        fprintf(repmat('\b',1,ll_collect))
        fprintf(repmat('\b',1,ll_loop1))
        
        clear f_temp plot_vars_temp data_temp
        
    end
    
    % Run through all remaining values
    if logical(chunk_rems)
        
        % Preallocate temp vars
        clear f_temp plot_vars_temp data_temp
        f_temp = zeros(chunk_rems,1);
        for i0 = 1:length(name)
            plot_vars_temp(chunk_rems).(name{i0}) = [];
        end
        data_temp(1:chunk_rems) = parallel.FevalFuture;

        % Display info
        ll_loop1 = fprintf('    Search: %1.0f-%1.0f / %1.0f\n',(chunk_whole*chunk_nums+1),(chunk_whole*chunk_nums+chunk_rems),num_iter_global);
        ll_loop2 = fprintf('    Search: %1.0f\n',(chunk_whole*chunk_nums+1));

        % Perform parfeval on remainder
        for chunk_rem = 1:chunk_rems

            % Search Display
            fprintf(repmat('\b',1,ll_loop2))
            ll_loop2 = fprintf('    Search: %1.0f\n',(chunk_whole*chunk_nums+chunk_rem));

            data_temp(chunk_rem) = parfeval(OPT.ppool,fh,2,x0((chunk_whole*chunk_nums+chunk_rem),:),BOD,CONST,OPT,VAR);
        end
        fprintf(repmat('\b',1,ll_loop2))

        % Collect results for chunk n
        ll_collect = fprintf('\n ---Collecting results from parpool---\n');
        for chunk_rem = 1:chunk_rems
            [index,cost,plt_vars] = fetchNext(data_temp);
            f_temp(index) = cost;
            plot_vars_temp(index) = plt_vars;
        end

        % Save the collected results into the total struct
        f( (chunk_whole*chunk_nums+1):num_iter_global ) = f_temp;
        plot_vars( (chunk_whole*chunk_nums+1):num_iter_global ) = plot_vars_temp;

        % Clear display
        fprintf(repmat('\b',1,ll_collect))
        fprintf(repmat('\b',1,ll_loop1))
        
    end

    % Clear display
    fprintf(repmat('\b',1,ll_RAM))
    fprintf(repmat('\b',1,ll_par))

else
    
    % Display info
    ll_loop = fprintf('    Search: %1.0f/%1.0f\n',1,OPT.MBH(count_MBH).N1_Outer);
    
    % Cost
    for i1 = 1:num_iter_global

        % Search Display
        fprintf(repmat('\b',1,ll_loop))
        ll_loop = fprintf('    Search: %1.0f/%1.0f\n',i1,num_iter_global); 

        % ***2.0 1st run of NLP***
        %Run the NLP problem solver to find x_star using the initial guess x0
        [f(i1),plot_vars(i1)] = feval(OPT.solver,x0(i1,:),BOD,CONST,OPT,VAR);

    end
    fprintf(repmat('\b',1,ll_loop))

end

neval = num_iter_global;



%% Check Feasibility

ll_feas = fprintf('\n ---Checking Feasibility---\n');
for i2 = 1:num_iter_global
    
    % ***1.0 Pre-allocate/Random point***
    f_best = [];
    x_star_best = [];
    N_not_improve = 0;
    
    % ***3.0 Check x_star feasibility***
    % ***MGALT NOTE***
    % Due to the multiple transfers and subsequently having to use an indirect
    % and direct method for solving, the variable x* is not going to be defined
    % as simplistic as it is within the psudo algorithm
    [feasible,~,~] = MGALT_MBH_isFeasible(BOD,CONST,OPT,VAR,plot_vars(i2),per_feas);
    
    % If there was feasibility
    if any(feasible)
%     if all(feasible)

        % A basin has been found, start the iteration process to explore it
        ll_MBH = fprintf('%s\n',disp_info);

        f_current = f(i2);
        x_current = x0(i2,:);

        f_best = f(i2);
        x_star_best = x0(i2,:);
        is_viable(i2) = true;
    else
        % No basin was found, exit the function to then try the next member
        f_best_return(i2) = f(i2);
        x_star_best_return(i2,:) = x0(i2,:);
        is_viable(i2) = false;
        continue
    end

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
        [feasible,~,~] = MGALT_MBH_isFeasible(BOD,CONST,OPT,VAR,plot_vars_prime,per_feas);
        if any(feasible)
%         if all(feasible)

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

    % Remove "Potential Basin Found"
    fprintf(repmat('\b',1,ll_MBH))


    % ***5.0 Find the problem solution***
    % Get the index of the best solution
    [~,index] = min(f_best);

    % Use the index to pull out the values
    f_best_return(i2) = f_best(index);
    x_star_best_return(i2,:) = x_star_best(index,1:end);

end
fprintf(repmat('\b',1,ll_feas))



end


