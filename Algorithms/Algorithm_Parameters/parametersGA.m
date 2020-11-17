function [options] = parametersGA(selection)
% FORM: [options] = parametersGA(selection)
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function call outputs the necessary struct object to fully
% |     define a Genetic Algorithm Island. A switch/case allows the user 
% |     to add in commonly used cases to have it as an easy return, and 
% |     more cases can easily be added. The function also allows for 
% |     custom inputs, which is then exported to be in the correct form 
% |     for the Island struct object.
% |
% |     -If the user desires another GA Island, all that is another call 
% |     to the "parametersGA" function. For example: the structure for 
% |     the settings of the 2nd island using Genetic Algorithm would be 
% |     named 
% |     "OPT.GA(1,2) = parametersGA('choice_here');"
% |
% |     -A set of parameters is located below in the section titled 'MISC'
% |     which describes the object struct and the different parts of the 
% |     object and what they do. For more details on the methods and 
% |     inputs here the user should refer to Appendix A of the 
% |     accompanying thesis.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -selection          (1,n)       [string]        [unitless]
% |         A used defined string which is used in a switch/case format
% |         to return predefined options
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -options            (1,n)       [assorted]      [unitless]
% |         The return struct with all the associated parameters
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     The parameters are defined as follows:
% |     -Npop:              (1,1)       [int]           [unitless]
% |         Number of members (solutions) in each population.
% |     -Ngen:              (1,1)       [int]           [unitless]
% |         Number of generations that will be evaluated.
% |     -pc:                (1,1)       [float]         [percent]
% |         Percent probability for a crossover to occur on a 0-1 scale.
% |     -pm:                (1,1)       [float]         [percent]
% |         Percent probablity for mutation to occur on a 0-1 scale.
% |     -gen_method:        (1,1)       [string]        [unitless]
% |         Generation method. This is how the algorithm will select 
% |         which members move to the next generation. Options are: 
% |         'total_random_replacement', 'tournament', 'natural_selection', 
% |         'thresholding', and 'weighted_random'. All options are case 
% |         sensitive.
% |     -N_keep:            (1,1)       [int]           [unitless]
% |         Number of solutions to keep between generations if the 
% |         generation method is 'natural_selection','thresholding', 
% |         or 'weighted_random'.
% |     -T:                 (1,1)       [int]           [unitless]
% |         Number of members to participate in each tournament if the 
% |         generation method is 'tournament'. MUST BE 2 OR GREATER.
% |     -elite:             (1,1)       [int]           [unitless]
% |         Number of elite solutions that automatically survive to 
% |         become part of the next generation. This is to prevent the 
% |         best solution from being lost so its recommended to be at 
% |         least 1.
% |     -threshold:         (1,1)       [int]           [uniltess]
% |         The cost threshold if the gereration method is 'thresholding'.
% |     -weight:            (1,1)       [string]        [unitless]
% |         How to calculate probabilities if the generation method is 
% |         'weighted_random'. Options are 'cost' and 'rank'.
% |     -mate_method:       (1,1)       [string]        [unitless]
% |         The mating method for two members. Options are 
% |         'uniform_crossover', 'random_crossover', and 'blending'.
% |     -cross_points:      (1,1)       [int]           [unitless]
% |         How many points will each member be allowed to be crossed 
% |         over with its mate. if the mate method is 'random_crossover'.
% |     -OB:                (1,1)       [float]         [percent]
% |         Out of bounds limit for blending (recommended < 20%).
% |
% |-----------------------------------------------------------------------



%% GA Parameters Selection

switch selection
    
    case {'75_30'}      % 75 popn and 30 gen
        
        options = struct(...
        'Npop',75,...                       % Number of members in each population
        'Ngen',30,...                       % Number of generations to evaluate
        'pc',0.80,...                       % Percent probability for a crossover to occur
        'pm',0.005,...                      % Percent probability for a mutation to occur
        'gen_method','tournament',...       % The generation method
        'N_keep',50,...                     % Number of solutions to keep between each generation
        'T',3,...                           % Number of members to participate in each tourneyment
        'elite',1,...                       % Number of elite solutions to survive
        'threshold',1000,...                % Cost threshold if 'gen_method' == 'threshold'
        'weight','cost',...                 % How to calcuate probabilities if 'gen_method' == 'weighted_random'
        'mate_method','random_crossover',...% The mating method for two members
        'cross_points',20,...               % How many points for the member to cross with the mate
        'OB',0.1);                          % Out of bounds limit for blending
        
    case {'100_30'}     % 100 popn and 30 gen

        options = struct(...
        'Npop',100,...                    	% Number of members in each population
        'Ngen',30,...                       % Number of generations to evaluate
        'pc',0.80,...                       % Percent probability for a crossover to occur
        'pm',0.005,...                      % Percent probability for a mutation to occur
        'gen_method','tournament',...       % The generation method
        'N_keep',50,...                     % Number of solutions to keep between each generation
        'T',3,...                           % Number of members to participate in each tourneyment
        'elite',1,...                       % Number of elite solutions to survive
        'threshold',1000,...                % Cost threshold if 'gen_method' == 'threshold'
        'weight','cost',...                 % How to calcuate probabilities if 'gen_method' == 'weighted_random'
        'mate_method','random_crossover',...% The mating method for two members
        'cross_points',20,...               % How many points for the member to cross with the mate
        'OB',0.10);                         % Out of bounds limit for blending

    case {'350_75'}     % 350 popn and 75 gen

        options = struct(...
        'Npop',350,...                     	% Number of members in each population
        'Ngen',75,...                       % Number of generations to evaluate
        'pc',0.80,...                       % Percent probability for a crossover to occur
        'pm',0.005,...                      % Percent probability for a mutation to occur
        'gen_method','tournament',...       % The generation method
        'N_keep',50,...                     % Number of solutions to keep between each generation
        'T',25,...                        	% Number of members to participate in each tourneyment
        'elite',1,...                       % Number of elite solutions to survive
        'threshold',1000,...                % Cost threshold if 'gen_method' == 'threshold'
        'weight','cost',...                 % How to calcuate probabilities if 'gen_method' == 'weighted_random'
        'mate_method','random_crossover',...% The mating method for two members
        'cross_points',20,...               % How many points for the member to cross with the mate
        'OB',0.10);                         % Out of bounds limit for blending
    
    otherwise           % Error

        errorPathDisplay();
        fprintf(2,'Incorrect GA selection.\n')
        return

end
       

end


