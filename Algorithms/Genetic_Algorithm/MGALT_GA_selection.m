function [members,fit] = MGALT_GA_selection(popn,f,OPT_algo)
% FORM: [members,fit] = MGALT_GA_selection(popn,f,opt_GA)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function selects the new members for the Genetic Algorithm
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -popn               (Npop,Nvar)	[float]         [unitless]
% |         Current members of the population
% |     -f                  (Npop,1)  	[float]         [unitless]
% |         The respective cost of each member in 'popn'
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         GA option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersGA.m"
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -members        	(Nkeep,Nvar)[float]         [unitless]
% |         Members selected to help produce the next generation
% |     -fit                (Nkeep,1)   [float]         [unitless]
% |         The respective cost of each member in 'members'
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup
[Npop,Nvar] = size(popn);

% Sensible Names to Variables
elite = OPT_algo.elite;          % Number of elite solutions to automatically enter next generation
N_keep = OPT_algo.N_keep;        % Number to select for mating
T = OPT_algo.T;                  % Number of members in each tournament
method = OPT_algo.gen_method;    % Selection method
threshold = OPT_algo.threshold;  % Threshold cost

% Preallocate
members = [];
fit = [];

% Sort Population from Best --> Worst
sorted_popn = popn*0;   % Pre-Allocation
[sorted_f,I] = sort(f);
for v = 1:Nvar
  sorted_popn(:,v) = popn(I,v);
end



%% Selection

switch method
    
    case {'natural_selection'}      % Natural Selection
        
        % Elite Members
        if elite
            members(1:elite,:) = sorted_popn(1:elite,:);
            fit(1:elite) = sorted_f(1:elite);
        end

        % Best N Members
        members(elite+1:N_keep+elite,:) = sorted_popn(elite+1:N_keep+elite,:);
        fit(elite+1:N_keep+elite) = sorted_f(elite+1:N_keep+elite);
        
    case {'tournament'}             % Tournament
        
    	% Which solutions have not been picked yet
        available = 1:length(sorted_f);
   
        % Keep "Elite" Members
        if elite
            for i = 1:elite
                members(i,:) = sorted_popn(i,:);
                fit(i) = sorted_f(i);
                available(i) = [];
            end
        end
   
        % Perform Tournament with Remaining members
        sel = elite + 1;
        while available
            % Select Gladiators
            if length(available) < T
                competitors = length(available);
            else
                competitors = T;
            end
            
            gladiator = zeros(competitors,Nvar);
            gladiator_f = zeros(1,competitors);
            
            for i = 1:competitors
                ind = randomNum(1,length(available),'int');
                gladiator(i,:) = sorted_popn(available(ind),:);
                gladiator_f(i) = sorted_f(available(ind));
                available(ind) = [];
            end
       
            % Fight to the Death
            [~,ind_victor] = min(gladiator_f);
            victor = gladiator(ind_victor,:);
            victor_f = gladiator_f(ind_victor);
       
            % Assign Vicotr to Selected Population
            members(sel,:) = victor;
            fit(sel) = victor_f;
            sel = sel + 1;
        end
        
    case {'thresholding'}           % Thresholding
        
    	% Keep "Elite" Members
        if elite
            for i = 1:elite
                members(i,:) = sorted_popn(i,:);
                fit(i) = sorted_f(i);
            end
        end
   
        % Find Solutions that Meet the Threshold
        i = elite + 1;
        while sorted_f(i) <= threshold
            members(i,:) = sorted_popn(i,:);
            fit(i) = sorted_f(i);
            i = i + 1;
            if i > Npop
                break
            end
        end
        
    case {'weighted_random'}        % Roulette Wheel
        
     	% Keep "Elite" Members
        if elite
            for i = 1:elite
                members(i,:) = sorted_popn(i,:);
                fit(i) = sorted_f(i);
            end
        end
   
        % Pre-Allocate Probability Array
        P = zeros(Npop-elite,1);
    
        % Assign Probabilities for Ranked Weighted Random
        if strcmpi(OPT_algo.weight,'rank')
            den = sum(1:Npop-elite);
            for member = 1:Npop-elite
                P(member) = (Npop-elite-member+1)/den;
            end

        % Assign Probabilities for Cost Weighted Random
        elseif strcmpi(OPT_algo.weight,'cost')
            C = -sorted_f(elite+1:end) + max(sorted_f);
            den = sum(C);
            if den ~= 0
                for member = 1:Npop-elite
                    P(member) = abs(C(member)/den);
                end
            else
                P = 1/(Npop-elite);
            end
        end
   
        for pick = 1:N_keep
            % "Spin" The Wheel
            spin = randomNum( 0,sum(P),'dec' );

            % March Through Probability Matrix Until The "Spun" Slice Is Hit
            for i = 1:length(P)
                if spin <= sum(P(1:i))
                    index = i;
                    members(pick+elite,:) = sorted_popn(i+elite,:); 
                    fit(pick+elite) = sorted_f(i+elite);
                    break
                end
            end
       
        P(index) = [];
        end
        
    otherwise
       
        errorPathDisplay();
        fprintf(2,'Incorrect generation method selected.\n')
        fprintf(2,'Check the documentation for the specific algorithm to see selection choices.\n\n')
        return
        
end



end


