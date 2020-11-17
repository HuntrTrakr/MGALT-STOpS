function [options] = parametersMBH(selection)
% FORM: [options] = parametersMBH(selection)
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function call outputs the necessary struct object to fully
% |     define a Monotonic Basin Hopping Island. A switch/case allows the 
% |     user to add in commonly used cases to have it as an easy return, 
% |     and more cases can easily be added. The function also allows for 
% |     custom inputs, which is then exported to be in the correct form 
% |     for the Island struct object.
% |
% |     -If the user desires another MBH Island, all that is another call 
% |     to the "parametersMBH" function. For example: the structure for 
% |     the settings of the 2nd island using Monotonic Basin Hopping 
% |     would be named 
% |     "OPT.MBH(1,2) = parametersMBH('choice_here');"
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
% |     -N1_Outer:          (1,1)       [int]           [unitless]
% |         Number of iterations to run for the first loop. These 
% |         iterations are the random population to explore the global 
% |         search space.
% |     -N1_Inner:          (1,1)       [int]           [unitless]
% |         Number of iterations to run within the basin finder in the 
% |         first loop. These iterations are only run if one of the 
% |         random populations from above are determined to be feasible 
% |         and hold a solution.
% |     -N2_Outer:          (1,1)       [int]           [unitless]
% |         Number of iterations to run for the second loop. These 
% |         iterations are the random population to explore the global 
% |         search space.
% |     -N2_Inner:          (1,1)       [int]           [unitless]
% |         Number of iterations to run within the basin finder in the 
% |         second loop. These iterations are only run if one of the 
% |         random populations from above are determined to be feasible 
% |         and hold a solution.
% |     -maxclst:           (1,1)       [int]           [unitless]
% |         When basins are found and sorted into clusters, how many 
% |         population arrays at a maximum shall be used within the 
% |         secondary solver loop. If the number "maxclst" is larger 
% |         than the amount of clusters found, it will decrease for 
% |         that iteration only.
% |     -per_feas:          (1,mig-1)	[float]         [unitless]
% |         A 1xn vector containing the percentage value for the 
% |         feasibility check. If, say a 500% and 200% value are 
% |         desired, the input is [5.00,2.00,...]. This percent value is 
% |         applied to the desired body's sphere of influence (SOI) value 
% |         in order to determine if the trajectory is going to terminate 
% |         around a desired region.
% |     -per_rand:          (1,mig-1)	[float]         [unitless]
% |         A 1xn vector containing the percentage value for the 
% |         perturbation of the initial population. If, say a 10% and 5% 
% |         value are desired, the input is [0.10,0.05,...]. This percent 
% |         value is applied to the members of the population to randomly 
% |         perturb them for the MBH algorithm scheme. Say one of the 
% |         initial population values is 250 and there is a 10% random 
% |         perturbation value applied to that value, a random number 
% |         will be generated between 225 and 275.
% |
% |-----------------------------------------------------------------------



%% MBH Parameters Selection

switch selection
    
    case {'3M-750_350'}     % 3 Mig: 750 loops and 350 iter

        options = struct(...
        'N1_Outer',750,...                  % Number of iter for the 1st run
        'N1_Inner',350,...                  % Number of iter per loop for the 1st
        'N2_Outer',200,...                  % Number of iter for the 2nd run
        'N2_Inner',50,...                   % Number of iter per loop for the 2nd
        'maxclst',30,...                    % Max number of clusters to find
        'per_feas',[5.00,2.50,1.25,100],...	% SOI to declare result feasible
        'per_rand',[0.10,0.05,0.02,0.01]); 	% Perturbation for MBH

    case {'3M-1200_500'}    % 3 Mig: 1200 loops and 500 iter

        options = struct(...
        'N1_Outer',1200,...                 % Number of iter for the 1st run
        'N1_Inner',500,...                  % Number of iter per loop for the 1st
        'N2_Outer',600,...                  % Number of iter for the 2nd run
        'N2_Inner',250,...                  % Number of iter per loop for the 2nd
        'maxclst',30,...                    % Max number of clusters to find
        'per_feas',[5.00,2.50,1.25,100],...	% SOI to declare result feasible
        'per_rand',[0.10,0.05,0.02,0.01]); 	% Perturbation for MBH

    case {'4M-500_200'}     % 4 Mig: 500 loops and 200 iter

        options = struct(...
        'N1_Outer',500,...                          % Number of iter for the 1st run
        'N1_Inner',200,...                          % Number of iter per loop for the 1st
        'N2_Outer',125,...                          % Number of iter for the 2nd run
        'N2_Inner',50,...                           % Number of iter per loop for the 2nd
        'maxclst',10,...                            % Max number of clusters to find
        'per_feas',[5.00,2.50,1.75,1.25,1.00],...	% SOI to declare result feasible
        'per_rand',[0.10,0.07,0.05,0.02,0.01]);     % Perturbation for MBH
    
    otherwise               % Error

        errorPathDisplay();
        fprintf(2,'Incorrect MBH selection.\n')
        return

end



end


