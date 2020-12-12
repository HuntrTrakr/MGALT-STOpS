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
% |     -f_best             (1,Nvar)    [float]         [unitless]
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

per_feas = OPT_algo.per_feas(num_mig);	% Feasibility percentage for SOI

search_index = [];
count_index = 1;

for i2 = 1:num_iter_global
    
    % ***3.0 Check x_star feasibility***
    [feasible,~,~] = MGALT_MBH_isFeasible(BOD,CONST,OPT,OPT_algo,VAR,plot_vars(i2),per_feas);
    
    if OPT_algo.feas_check(feasible)
        
        search_index(count_index) = i2;
        count_index = count_index+1;
        
        is_viable(i2) = true;

    else
        
        % No basin was found, try the next member
        f_best_return(i2) = f(i2);
        x_star_best_return(i2,:) = x0(i2,:);
        is_viable(i2) = false;
        
    end
    
end

% If nothing feasible was found, return
if size(search_index,2) == 0
    return
end

fprintf(repmat('\b',1,ll_feas))



%% If Feasibility, Check Basins

if loop == 1
    disp_info = '\n ---Potential Basin Found---';
elseif loop == 2
    disp_info = '\n ---Exploring Minima---';
end

% A basin has been found, start the iteration process to explore it
ll_MBH = fprintf('%s\n',disp_info);

if OPT.parallel
    
    % Preallocate
    fh = @MGALT_MBH_basin;
    b_data(1:size(search_index,2)) = parallel.FevalFuture;
    
    % Display info
    ll_par = fprintf('\n    Parallel Processing   \n');
    ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n');
    
    % feval
    for i3 = 1:size(search_index,2)
        b_data(i3) = parfeval(OPT.ppool,fh,2,BOD,CONST,OPT,OPT_algo,VAR,...
            f(search_index(i3)),x0(search_index(i3),:),num_mig,count_MBH,loop);
    end
    
     % Collect results as they become available
    for i4 = 1:size(search_index,2)
        [p_index,p_f,p_x] = fetchNext(b_data);
        p_f_best(p_index) = p_f;
        p_x_best(p_index,:) = p_x;
    end
    
    % Put the data into the correct index for func return
    for i5 = 1:size(search_index,2)
        f_best_return(search_index(i5)) = p_f_best(i5);
        x_star_best_return(search_index(i5),:) = p_x_best(i5,:);
    end

    % Clear display
    fprintf(repmat('\b',1,ll_RAM))
    fprintf(repmat('\b',1,ll_par))
    
else

    for i3 = 1:size(search_index,2)

        [f_best_return(search_index(i3)),x_star_best_return(search_index(i3),:)] = ...
            MGALT_MBH_basin(BOD,CONST,OPT,OPT_algo,VAR,...
            f(search_index(i3)),x0(search_index(i3),:),num_mig,count_MBH,loop);

    end

end

% Remove "Potential Basin Found"
fprintf(repmat('\b',1,ll_MBH))



end


