function [opt_island] = optionsIsland(selection)
% FORM: [opt_island] = optionsIsland(selection)
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function call outputs the necessary struct object to fully
% |     define the Island and algorithm connections. A switch/case allows 
% |     the user to add in commonly used cases to have it as an easy 
% |     return, and more cases can easily be added. The function also 
% |     allows for custom inputs, which is then exported to be in the 
% |     correct form for the Island struct object.
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
% |     -opt_Isl            (1,n)       [assorted]      [unitless]
% |         The return object with all the associated parameters
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     The parameters are defined as follows:
% |     -Nmig               (1,1)       [int]           [unitless]
% |         How many migrations
% |     -isl_list           (n,1)       [cell]          [string]
% |         Which algorithms to use and in which order
% |         Ex: {'MBH'}     Ex: {'GA';'PSO'}
% |     -isl_conn           (n,n)       [int]           [unitless]
% |         The connection matrix between the islands
% |         Ex: [1]         Ex: [1,1; 1,1]
% |     -rep_pol           (n,1)       [cell]          [string]
% |         The replacement policy for the islands
% |         Ex: {'best_n'}  Ex: {'best_n';'best_n'}
% |     -rep_opt           (n,1)       [int]           [unitless]
% |         Replacement options for each island
% |         Ex: [10]        Ex: [10,10]
% |     -sel_pol           (n,1)       [cell]          [string]
% |         Selection policy for each island
% |         Ex: {'natural_selection'}
% |         Ex: {'natural_selection';'natural_selection'}
% |     -sel_opt           (1,n)       [int]           [unitless]
% |         Number of solutions to select for each selection policy
% |         Ex: [10]        Ex: [10,10]
% |
% |-----------------------------------------------------------------------



%% Selection

opt_island = struct;

switch selection
    
    case {'2M-all'}
        
        opt_island.Nmig = 2;                           	% How many migrations
        opt_island.isl_list = {'GA';'DE';'PSO';'MBH'};	% Which algorithms to use and in which order
        opt_island.isl_conn = [1,1,1,1; 1,1,1,1;...
            1,1,1,1; 1,1,1,1];                          % The connection matrix between the islands
        opt_island.rep_pol = {'best_n';'best_n';...
            'best_n';'best_n'};                         % The replacement policy for the islands
        opt_island.rep_opt = [10,10,10,10];          	% Replacement options for each island
        opt_island.sel_pol = {'natural_selection';...
            'natural_selection';...
            'natural_selection';...
            'natural_selection'};                       % Selection policy for each island
        opt_island.sel_opt = [10,10,10,10];           	% Number of solutions to select for each selection policy
        
    case {'3M-all'}
        
        opt_island.Nmig = 3;                           	% How many migrations
        opt_island.isl_list = {'GA';'DE';'PSO';'MBH'};	% Which algorithms to use and in which order
        opt_island.isl_conn = [1,1,1,1; 1,1,1,1;...
            1,1,1,1; 1,1,1,1];                          % The connection matrix between the islands
        opt_island.rep_pol = {'best_n';'best_n';...
            'best_n';'best_n'};                         % The replacement policy for the islands
        opt_island.rep_opt = [10,10,10,10];          	% Replacement options for each island
        opt_island.sel_pol = {'natural_selection';...
            'natural_selection';...
            'natural_selection';...
            'natural_selection'};                       % Selection policy for each island
        opt_island.sel_opt = [10,10,10,10];           	% Number of solutions to select for each selection policy
        
    case {'4M-all'}
        
        opt_island.Nmig = 4;                           	% How many migrations
        opt_island.isl_list = {'GA';'DE';'PSO';'MBH'};	% Which algorithms to use and in which order
        opt_island.isl_conn = [1,1,1,1; 1,1,1,1;...
            1,1,1,1; 1,1,1,1];                          % The connection matrix between the islands
        opt_island.rep_pol = {'best_n';'best_n';...
            'best_n';'best_n'};                         % The replacement policy for the islands
        opt_island.rep_opt = [10,10,10,10];          	% Replacement options for each island
        opt_island.sel_pol = {'natural_selection';...
            'natural_selection';...
            'natural_selection';...
            'natural_selection'};                       % Selection policy for each island
        opt_island.sel_opt = [10,10,10,10];           	% Number of solutions to select for each selection policy
        
    case {'5M-all'}
    
        opt_island.Nmig = 5;                           	% How many migrations
        opt_island.isl_list = {'GA';'DE';'PSO';'MBH'};	% Which algorithms to use and in which order
        opt_island.isl_conn = [1,1,1,1; 1,1,1,1;...
            1,1,1,1; 1,1,1,1];                          % The connection matrix between the islands
        opt_island.rep_pol = {'best_n';'best_n';...
            'best_n';'best_n'};                         % The replacement policy for the islands
        opt_island.rep_opt = [10,10,10,10];          	% Replacement options for each island
        opt_island.sel_pol = {'natural_selection';...
            'natural_selection';...
            'natural_selection';...
            'natural_selection'};                       % Selection policy for each island
        opt_island.sel_opt = [10,10,10,10];           	% Number of solutions to select for each selection policy
    
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect algorithm input for Island option.\n')
        return
        
end


%The number of islands
opt_island.Nisl = length(opt_island.isl_list);



end


