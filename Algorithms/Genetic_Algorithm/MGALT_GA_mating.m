function [child1,child2] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2)
% FORM: [child1,child2] = MGALT_GA_mating(BOD,OPT,OPT_algo,VAR,parent1,parent2)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function allows two different parents to mate and produce 
% |     two different offspring. Used in the Genetic Algorithm (GA)
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
% |     -parent1          	(1,Nvar)	[float]         [unitless]
% |         First parent member of the population
% |     -parent2            (1,Nvar)  	[float]         [unitless]
% |         Second parent member of the population
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -child1             (1,Nvar)	[float]         [unitless]
% |         First child member of the new population
% |     -child2             (1,Nvar)	[float]         [unitless]
% |         Second child member of the new population
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

Nvar = length(parent1);
child1 = parent1;
child2 = parent2; % In case neither crossover nor mutation occurs



%% Mating

% Sensible Names to Variables
high = VAR.high;
low = VAR.low;
bin = VAR.bin;
pc = OPT_algo.pc;
pm = OPT_algo.pm;

switch OPT_algo.mate_method
    
    case {'random_crossover','uniform_crossover'}
        
        % For 3 digit availability, slide variables over 3 digits
        fac = 1000;
        p_bin1 = [];
        p_bin2 = [];
   
        % If there is a negative number, slide all variables up by that number
        min_num = min([parent1,parent2]);
        if min_num < 0
            parent1 = parent1 + abs(min_num);
           	parent2 = parent2 + abs(min_num);
        end

        % Get all binary arrays; Make matching arrays equal in length
        for i = 1:Nvar
            p1 = dec2bin(parent1(i)*fac);
            
            if isinf(parent1(i)*fac) || isinf(parent2(i)*fac)
                disp('u')
            end
            
            p2 = dec2bin(parent2(i)*fac);
            
            while length(p1) < length(p2)
                p1 = [ '0', p1];
            end
            while length(p2) < length(p1)
                p2 = [ '0', p2];
            end
            
            l(i) = length(p1); %#ok<AGROW>
            p_bin1 = [ p_bin1, p1 ];
            p_bin2 = [ p_bin2, p2 ]; %#ok<AGROW>
        end
        kid1 = p_bin1;
        kid2 = p_bin2;

        % Perform the bit crossover
        N = length(p_bin1);
        if OPT_algo.cross_points > (N-1)
            OPT_algo.cross_points = N-1; 
        end
        if strcmp(OPT_algo.mate_method,'uniform_crossover')
            OPT_algo.cross_points = N-1; 
        end
        available = 1:1:(N-1);
        
        for i = 1:OPT_algo.cross_points % Determine crossover cutoff points
            choice = randomNum(1,length(available),'int');
            cut(i) = available(choice); %#ok<AGROW>
            available(choice) = [];
        end
        
        cut(end+1) = 0;
        cut = sort(cut);
   
        for i = 1:length(cut)-1
            array = (cut(i)+1):cut(i+1);
            if rand <= pc
                kid1(array) = p_bin2(array); 
            end
            if rand <= pc
                kid2(array) = p_bin1(array); 
            end
        end
   
        % Turn the binary strings back into variable arrays
        for i = 1:Nvar
            array = 1:l(i);
            if i ~= 1
                array = array + sum(l(1:i-1));
            end
            child1(i) = bin2dec(kid1(array))/fac;
            child2(i) = bin2dec(kid2(array))/fac;
        end
   
        % If the number were slid to avoid negatives, slide them back
        if min_num < 0
            child1 = child1 - abs(min_num);
            child2 = child2 - abs(min_num);
        end
        
    case {'blending'}
        
    	child1 = parent1; child2 = parent2;
        if parent1 ~= parent2 
            for i = 1:Nvar
                if rand <= pc
                    beta1 = randomNum(-OPT_algo.OB,1+OPT_algo.OB,'dec');
                   	child1(i) = beta1*parent1(i) + (1-beta1)*parent2(i);
                end
                if rand <= pc
                    beta2 = randomNum(-OPT_algo.OB,1+OPT_algo.OB,'dec');
                   	child2(i) = beta2*parent1(i) + (1-beta2)*parent2(i);
                end
            end
   
        % If both parents are the same, there is a chance to 'bump' the
        % solution a little bit to see if it improves, instead of just
        % having both children be the same as the parents
        else
            range = high - low;
            for i = 1:Nvar
                if rand <= pc
                    beta1 = randomNum(-OPT_algo.OB,OPT_algo.OB,'dec');
                  	child1(i) = (beta1*range(i)) + parent1(i);
                end
                if rand <= pc
                    beta2 = randomNum(1-OPT_algo.OB,1+OPT_algo.OB,'dec');
                  	child2(i) = (beta2*range(i)) + parent1(i);
                end
            end
        end
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect mating method selected.\n')
        fprintf(2,'Check the documentation for the specific algorithm to see selection choices.\n\n')
        return
        
end



%% Mutation

% Mutate All Variables
if rand <= pm
   for i = 1:Nvar
      child1(i) = (high(i)-low(i))*rand + low(i);
   end
end
if rand <= pm
   for i = 1:Nvar
      child2(i) = (high(i)-low(i))*rand + low(i);
   end
end



%% Ensure Bounds Not Broken and Binary Stays Binary

for i = 1:Nvar
   if child1(i) > high(i)
       child1(i) = high(i);
   end
   
   if child1(i) < low(i)
       child1(i) = low(i);
   end
   
   if child2(i) > high(i)
       child2(i) = high(i); 
   end
   
   if child2(i) < low(i)
       child2(i) = low(i);
   end
   
   if bin(i)
       child1(i) = round(child1(i));
       child2(i) = round(child2(i));
   end
end



%% Correct for ToF 
% Added by Malloy to account for the incorrect ToF variables between the
% members

child1 = MGALT_fixToF(BOD,OPT,VAR,child1);
child2 = MGALT_fixToF(BOD,OPT,VAR,child2);



end


