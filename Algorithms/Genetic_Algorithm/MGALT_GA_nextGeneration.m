function [new_popn] = MGALT_GA_nextGeneration(BOD,OPT,OPT_algo,VAR,old_popn,f)
% FORM: [new_popn] = MGALT_GA_nextGeneration(BOD,OPT,OPT_algo,VAR,old_popn,f)
%
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
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         GA option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersGA.m"
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -old_popn         	(Npop,Nvar)	[float]         [unitless]
% |         The old population members which will be used to create the 
% |         new population members
% |     -f                  (Npop,1)  	[float]         [unitless]
% |         The respective cost of each member in 'old_popn'
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -new_popn           (Npop,Nvar)	[float]         [unitless]
% |         The new population members which generated from the old 
% |         population members and random crossovers
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

% Assign Sensible Names
[Npop,Nvar] = size(old_popn);
new_popn = old_popn*0; % Pre-Allocation
elite = OPT_algo.elite;



%% Replacement

switch OPT_algo.gen_method
    
    case {'total_random_replacement'}   % Total Random Replacement
        
        for member = 1:2:Npop
        	% Select 2 Parents at Random
            parent1 = old_popn(randomNum(1,Npop,'int'),:);
            parent2 = old_popn(randomNum(1,Npop,'int'),:);

            % Generate 2 Offspring
            [new_popn(member,:),kid] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2);
            if (member+1) <= Npop
                new_popn(member+1,:) = kid;
            end
        end
 
    case {'natural_selection'}          % Natural Selection
        
        % Select Mating population
        [sel_popn,sel_f]= MGALT_GA_selection(old_popn,f,OPT_algo);
   
        % Only Keep Elite Members
        if elite
            new_popn(1:elite,:) = sel_popn(1:elite,:);
        end
   
        % Fill Rest of Population
        for member = (elite+1):2:Npop
            % Select 2 "Elite" Parents at Random
            parent1 = sel_popn(randomNum(1,length(sel_f),'int'),:);
            parent2 = sel_popn(randomNum(1,length(sel_f),'int'),:);
      
            % Generate 2 Offspring
            [new_popn(member,:),kid] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2);
            if (member+1) <= Npop
                new_popn(member+1,:) = kid; 
            end
        end
        
    case {'tournament'}                 % Tournament
        
        % Select Mating population
        [sel_popn,sel_f] = MGALT_GA_selection(old_popn,f,OPT_algo);
   
        % Only Keep Elite Members
        if elite
            new_popn(1:elite,:) = sel_popn(1:elite,:);
        end
   
        % Fill Rest of Population
        for member = (elite+1):2:Npop
            % Select 2 "Elite" Parents at Random
            parent1 = sel_popn(randomNum(1,length(sel_f),'int'),:);
            parent2 = sel_popn(randomNum(1,length(sel_f),'int'),:);
      
            % Generate 2 Offspring
            [new_popn(member,:),kid] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2);
            if (member+1) <= Npop
                new_popn(member+1,:) = kid; 
            end
        end
        
    case {'thresholding'}               % Thresholding
        
        % Select Mating population
        [sel_popn,sel_f] = MGALT_GA_selection(old_popn,f,OPT_algo);

        % Only Keep Elite Members
        if elite
            new_popn(1:elite,:) = sel_popn(1:elite,:);
        end

        % Re-Initialize Random Population If Criteria Never Met
        if isempty(sel_f) || length(sel_f) == elite
            new_popn(elite+1:Npop,1:Nvar) = rand(Npop-elite,Nvar);
            for p = elite+1:Npop
                for v = 1:Nvar
                    if ~VAR.bin(v) % non-binary variables
                        new_popn(p,v) = (VAR.high(v)-VAR.low(v))*new_popn(p,v) + VAR.low(v);
                    else % binary variables
                        new_popn(p,v) = round(new_popn(p,v));
                    end
                end
            end
        end

        % Fill Rest of Population
        if length(sel_f) > elite 
            for member = (elite+1):2:Npop
                % Select 2 "Elite" Parents at Random
                parent1 = sel_popn(randomNum(1,length(sel_f),'int'),:);
                parent2 = sel_popn(randomNum(1,length(sel_f),'int'),:);

                % Generate 2 Offspring
                [new_popn(member,:),kid] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2);
                if (member+1) <= Npop
                    new_popn(member+1,:) = kid; 
                end
            end
        end
        
    case {'weighted_random'}            % Weighted Random
        
        % Select Mating population
        [sel_popn,sel_f] = MGALT_GA_selection(old_popn,f,OPT_algo);
   
        % Only Keep Elite Members
        if elite
            new_popn(1:elite,:) = sel_popn(1:elite,:);
        end
   
        % Fill Rest of Population
        for member = (elite+1):2:Npop
            % Select 2 "Elite" Parents at Random
            parent1 = sel_popn(randomNum(1,length(sel_f),'int'),:);
            parent2 = sel_popn(randomNum(1,length(sel_f),'int'),:);
      
            % Generate 2 Offspring
            [new_popn(member,:),kid] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2);
            if (member+1) <= Npop
                new_popn(member+1,:) = kid; 
            end
        end
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect replacement method selected.\n')
        fprintf(2,'Check the documentation for the specific algorithm to see selection choices.\n\n')
        return
        
end



end


