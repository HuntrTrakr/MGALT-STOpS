function [new_popn,new_f,nfeval] = MGALT_DE_nextGeneration(BOD,CONST,OPT,...
    OPT_algo,VAR,old_popn,old_f,nfeval,count_DE)
% FORM: [new_popn,new_f,nfeval] = MGALT_DE_nextGeneration(BOD,CONST,OPT,...
%       OPT_algo,VAR,old_popn,old_f,nfeval,count_DE)
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function creates the next generation of DE popn from the 
% |     old member population and slight perturbations to the old members. 
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
% |         DE option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersDE.m"
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -old_popn         	(Npop,Nvar)	[float]         [unitless]
% |         The old population members which will be used to create the 
% |         new population members
% |     -old_f              (Npop,1)  	[float]         [unitless]
% |         The respective cost of each member in 'old_popn'
% |     -nfeval           	(1,1)       [int]           [unitless]
% |         The current number of iterations for this DE object
% |     -count_DE           (1,1)       [int]           [unitless]
% |         The current DE island number
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -new_popn           (Npop,Nvar)	[float]         [unitless]
% |         The new population members which generated from the old 
% |         population members and random perturbations
% |     -new_f              (Npop,1)    [float]         [unitless]
% |         The respective cost of each member in 'new_popn'
% |     -nfeval           	(1,1)       [int]           [unitless]
% |         The current number of iterations for this DE object
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

[Npop,Nvar] = size(old_popn);

% Pre-Allocation
trial_popn = old_popn;
new_popn = old_popn*0;
trial_f = old_f; 
new_f = old_f*0;



%% Base Vector Selection

base_choice = (1:Npop)*0;

switch OPT_algo.sel_method
    
    case {'random'}             % Random
        
    	% Select Base Vectors So That r0 ~= i
        available = 1:1:Npop;
        for i = 1:Npop
            choice = randomNum(1,length(available),'int');
            base_choice(i) = available(choice);
            available(choice) = [];
        end
        
        % If Any r0 == i, Fix It
        ctr = 0;
        
        while any(base_choice == 1:Npop)
            I = find(base_choice == 1:Npop);
            
            for i = 1:length(I)
                temp = base_choice(I(i));
                
                if I(i) > 1
                    base_choice(I(i)) = base_choice(I(i)-1);
                    base_choice(I(i)-1) = temp;
                else
                    base_choice(I(i)) = base_choice(I(i)+1);
                   	base_choice(I(i)+1) = temp;
                end
                
            end
            
            ctr = ctr + 1;
            
            if ctr > 11
                available = 1:1:Npop;
                
                for i = 1:Npop
                    choice = randomNum(1,length(available),'int');
                    base_choice(i) = available(choice);
                    available(choice) = [];
                end
                
            end
            
        end
        
        base_vec = old_popn(base_choice,:);
        
    case {'best_so_far'}        % Best so Far
        
        [~, ind] = min(old_f);
        base_choice = base_choice + ind;
        base_vec = old_popn(base_choice,:);
        
    case {'random_best_blend'}  %Random Best Blend
        
       base_vec = old_popn*0;
        [~, ind] = min(old_f);
        best_so_far = old_popn(ind,:);
        
        % Select Base Vectors So That r0 ~= i
        available = 1:1:Npop;
        for i = 1:Npop
            choice = randomNum(1,length(available),'int');
            base_choice(i) = available(choice);
            available(choice) = [];
        end
        
        % If Any r0 == i, Fix It
        ctr = 0;
        
        while any(base_choice == 1:Npop)
            I = find(base_choice == 1:Npop);
            
            for i = 1:length(I)
                temp = base_choice(I(i));
                
                if I(i) > 1
                    base_choice(I(i)) = base_choice(I(i)-1);
                   	base_choice(I(i)-1) = temp;
                else
                    base_choice(I(i)) = base_choice(I(i)+1);
                 	base_choice(I(i)+1) = temp;
                end
                
            end
            
            ctr = ctr + 1;
            
            if ctr > 11
                available = 1:1:Npop;
                
                for i = 1:Npop
                    choice = randomNum(1,length(available),'int');
                    base_choice(i) = available(choice);
                    available(choice) = [];
                end
                
            ctr = 0;
            
            end
            
        end
        
        random_choice = old_popn(base_choice,:);
   
        for i = 1:Npop
            base_vec(i,:) = random_choice(i,:) + rand*(best_so_far - random_choice(i,:));
        end 
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect selection method selected.\n')
        fprintf(2,'Check the documentation for the specific algorithm to see selection choices.\n\n')
        return
        
end



%% Generate All Trial Vectors

% Run Through All Target Vectors
if OPT.parallel

    % Display Info
    ll_par = fprintf('\n    Parallel Processing   \n');
    ll_RAM = fprintf(2,'    Watch CPU/RAM usage   \n\n');
    
    % Preallocate
    data(1:OPT.DE(count_DE).Npop) = parallel.FevalFuture;
    fh = str2func(string(OPT.solver));      % https://www.mathworks.com/help/matlab/ref/str2func.html
    
    for tv = 1:Npop
        % Generate Target Vectors
        [trial_popn(tv,:)] = MGALT_DE_Target(BOD,OPT,OPT_algo,VAR,...
            old_popn,base_choice,base_vec,Npop,Nvar,tv);
    end

    % Cost
    for tv = 1:Npop
        % Run trial solution
        data(tv) = parfeval(OPT.ppool,fh,1,trial_popn(tv,:),BOD,CONST,OPT,VAR);
    end
    
    % Collect results as they become available
    for tv = 1:Npop
        [index,cost] = fetchNext(data);
        trial_f(index) = cost;
    end
    
    fprintf(repmat('\b',1,ll_RAM))
    fprintf(repmat('\b',1,ll_par))
        

else

    for tv = 1:Npop

        % Generate Target Vectors
        [trial_popn(tv,:)] = MGALT_DE_Target(BOD,OPT,OPT_algo,VAR,...
            old_popn,base_choice,base_vec,Npop,Nvar,tv);

        % Run trial solution
        [trial_f(tv),~] = feval(OPT.solver,trial_popn(tv,:),BOD,CONST,OPT,VAR);

    end

end
    

nfeval = nfeval + Npop;
all_popn = [old_popn;trial_popn];
all_f = [old_f;trial_f];
 


%% Determine Which Members Survive

sorted_popn = all_popn*0; % Pre-Allocation

switch OPT_algo.surv_method
    
    case {'natural_selection'}  % Natural Selection
        
        [sorted_f,I] = sort(all_f);
        
        for v = 1:Nvar
            sorted_popn(:,v) = all_popn(I,v);
        end

        new_popn = sorted_popn(1:Npop,:);
        new_f= sorted_f(1:Npop);
        
    case {'tournament'}         % Tournament
        
        wins = all_f*0; % Pre-Allocation
        
        for i = 1:(2*Npop)
            available = [1:(i-1) , (i+1):(2*Npop)];
            tally = 0;
            
            for comp = 1:OPT_algo.T
                c = randomNum(1,length(available),'int');
                competitor = available(c);
                
                if all_f(competitor) <= all_f(i)
                    tally = tally + 1;
                end
                
                available(c) = [];
            end
            
            wins(i) = tally;
        end
        
        [~ , I] = sort(wins);
        
        for v = 1:Nvar
            sorted_popn(:,v) = all_popn(I,v);
        end
        
        sorted_f = all_f(I);

        new_popn = sorted_popn(1:Npop,:);
        new_f    = sorted_f   (1:Npop)  ;
        
    case {'weighted_random'}    % Weighted Random
        
        % Pre-Allocate Probability Array
        P = zeros(2*Npop,1);
        
        % Probability Assignment
        switch OPT_algo.weight
            
            case {'rank'}
                
                den = sum(1:(2*Npop));
            
                for member = 1:(2*Npop)
                    P(member) = ((2*Npop)-member+1)/den;
                end
                
            case {'cost'}
                
                normalizer = min(all_f);
            
                for member = 1:(2*Npop)
                    C(member) = all_f(member) - normalizer; %#ok<AGROW>
                end

                den = sum(C);

                if den ~= 0
                    for member = 1:(2*Npop)
                        P(member) = abs(C(member)/den);
                    end
                else
                    P = P + 1/(2*Npop);
                end
                
            otherwise
            
                errorPathDisplay();
                fprintf(2,'Incorrect difference vector selected.\n');
                disp('Incorrect difference vector selected.')
                return
            
        end
        
        available = (1:1:(2*Npop))';
        
        for i = 1:Npop
            
            % "Spin" The Wheel
            spin = randomNum( 0,sum(P),'dec' );
            
            % March Through Probability Matrix Until The "Spun" Slice Is Hit
            for j = 1:length(P)
                
                if spin <= sum(P(1:j))
                    new_popn(i,:) = all_popn(available(j),:); 
                    new_f   (i)   = all_f   (available(j))  ; 
                    P(j) = []; 
                    available(j) = [];
                    break;
                end
                
            end
            
        end
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect survival method selected.\n')
        return 
        
end



end


