function [trial_popn] = MGALT_DE_Target(BOD,OPT,OPT_algo,VAR,old_popn,...
    base_choice,base_vec,Npop,Nvar,tv)
% FORM: [trial_popn] = MGALT_DE_Target(BOD,OPT,OPT_algo,VAR,old_popn,...
%       base_choice,base_vec,Npop,Nvar,tv)
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function creates the target vector for DE_nextGeneration
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
% |     -base_choice     	(1,1)       [float]         [unitless]
% |         Index number for the base choice
% |     -base_vec           (1,1)       [float]         [unitless]
% |         Index number for the base vector selection
% |     -Npop            	(1,1)       [int]           [unitless]
% |         The size of the population
% |     -Nvar             	(1,1)       [int]           [unitless]
% |         The number of variables per member
% |     -tv                 (1,1)       [int]           [unitless]
% |         The current index in "Npop"
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -trial_popn     (1,old_popn)    [float]         [unitless]
% |         The trial population
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

pc = OPT_algo.pc;   	% Probability of Crossover
F = OPT_algo.F;      	% Scaling Factor

% Choose Target Vector
target_vec = old_popn(tv,:);



%% Members

% Choose Two Random Population Members So That r1 ~= r2 ~= i
choice = randomNum(1,Npop,'int');
while (choice == tv) || (choice == base_choice(tv))
    choice = randomNum(1,Npop,'int');
end
choice1 = choice;

choice = randomNum(1,Npop,'int');
while (choice == tv) || (choice == base_choice(tv)) || (choice == choice1)
    choice = randomNum(1,Npop,'int');
end
choice2 = choice;

% Compute Weighted Difference Vector
switch OPT_algo.F_method

    case {'constant'}

        difference_vec = F*(old_popn(choice1,:) - old_popn(choice2,:));

    case {'jitter'}

        for i = 1:Nvar
            difference_vec(i) = randomNum(F(1),F(2),'dec')*(old_popn(choice1,i) - old_popn(choice2,i));
        end

    case {'dither'}

        difference_vec = randomNum(F(1),F(2),'dec')*(old_popn(choice1,:) - old_popn(choice2,:));

    otherwise

        errorPathDisplay();
        fprintf(2,'Incorrect difference vector selected.\n')
        return

end

% Add to Base Vector
mutant_vec = base_vec(tv,:) + difference_vec;

% Crossover
trial_popn = target_vec;
for i = 1:Nvar
    if rand <= pc
        trial_popn(1,i) = mutant_vec(i);
    end

    % Ensure bounds not breached
    if trial_popn(1,i) > VAR.high(i)
        trial_popn(1,i) = VAR.high(i);
    elseif trial_popn(1,i) < VAR.low (i)
        trial_popn(1,i) = VAR.low(i);
    elseif VAR.bin(i)
        trial_popn(1,i) = round(trial_popn(1,i));
    end
end

% Ensure ToF is not borken
trial_popn(1,:) = MGALT_fixToF(BOD,OPT,VAR,trial_popn(1,:));
    


end


