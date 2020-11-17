function [popn,bee] = genPopn(algorithm,OPT,OPT_algo,VAR)
% FORM: [popn,bee] = genPopn(algorithm,OPT,opt_algo,VAR)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function designed to randomly generate an initial population for 
% |     one of the algorithms. Each algorithms has a specificway to 
% |     generate the population, which is specified in the switch/case 
% |     by the algorithm string. This funciton also ensures that the Tof 
% |     bounds are not broken from the start, so T1_JD+T1_ToF == T2_JD
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -algorithm          (1,n)       [string]        [unitless]
% |         The solver (cost function) handle which is used to determine 
% |         the correct solver to use. This string is used in an ode45 
% |         function call and is case specific for the function.m name
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         Algo option parameters
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -popn           (Nshared,Nvar)	[float]         [unitless]
% |         A randomely generated array of members from the initalization 
% |         of an island
% |     -bee                (Npop,1)	[struct]       	[unitless]
% |         Contains cost, speed position, ect of every bee in the 
% |         population for the new generation
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initials

Nvar = size(VAR.bin,2);       	% Number of Variables




%% Pick which solver/algorithm to generate for

switch OPT.solver
    
    case {'LT_DIR_FSM_2D','LT_IN_FSM_2D'}
        
        segment = Nvar/VAR.transfers;	% Number of segments for the population
        
        switch algorithm
    
            case {'GA','DE'}

                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.Npop,Nvar); %population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end

                % Break the popn into transfer segments
                % ie, popn is 75x44 with 2 transfers, so the first stansfer 
                % segment is popn(:,1:21) and the second segment is popn(:,22:end)
                % This needs to account for the time of flight, as the second
                % population segment can't start before the first one, so the
                % transfer departure JD must be adjusted
                for i3 = 1:OPT_algo.Npop
                    for i4 = 1:(VAR.transfers-1)
                        start = popn(i3,(segment*i4)-segment+1);
                        tof = popn(i3,segment*i4);
                        popn(i3,(segment*i4)+1) = start+tof;
                    end
                end

            case {'PSO'}

                % Initials
                Vmax = OPT_algo.vmax*(VAR.high-VAR.low);

                % Pre-Allocate
                popn = zeros(OPT_algo.Npop,Nvar);
                bee(1:OPT_algo.Npop) = struct('pos',zeros(1,Nvar), 'vel',zeros(1,Nvar), 'f',0,...
                    'f_p',999999999999999, 'p',zeros(1,Nvar), 'f_g',999999999999999, 'g',zeros(1,Nvar)); 

                % Randomize the Values
                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        bee(i1).pos(i2) = randomNum(VAR.low(i2),VAR.high(i2),'dec');
                        bee(i1).vel(i2) = randomNum(-Vmax(i2),Vmax(i2),'dec');
                        if VAR.bin(i2)
                            bee(i1).pos(i2) = round(bee(i1).pos(i2));
                        end
                    end
                end

                % Break the popn into transfer segments
                % ie, popn is 75x40 with 2 transfers, so the first stansfer 
                % segment is popn(:,1:20) and the second segment is popn(:,21:end)
                for i3 = 1:OPT_algo.Npop
                    for i4 = 1:(VAR.transfers-1)
                        start = bee(i3).pos((segment*i4)-segment+1);
                        tof = bee(i3).pos(segment*i4);
                        bee(i3).pos((segment*i4)+1) = start+tof;
                    end
                    popn(i3,:) = bee(i3).pos;
                end

            case {'MBH'}

                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.N1_Outer,Nvar); % Population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.N1_Outer
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2); 
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end

                % Break the popn into transfer segments
                % ie, popn is 75x40 with 2 transfers, so the first stansfer 
                % segment is popn(:,1:20) and the second segment is popn(:,21:end)
                for i3 = 1:OPT_algo.N1_Outer
                    for i4 = 1:(VAR.transfers-1)
                        start = popn(i3,(segment*i4)-segment+1);
                        tof = popn(i3,segment*i4);
                        popn(i3,(segment*i4)+1) = start+tof;
                    end
                end

            otherwise

                errorPathDisplay();
                errorAlgorithm();
                return

        end
        
    case {'MGALT_DIR_FBSM_2D'}
        
       	switch algorithm
            
            case {'GA','DE'}
                
                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.Npop,Nvar); %population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end
                
                % Break the popn into transfer segments
                if VAR.transfers > 1
                    
                    % How many segments is the transfers broken into
                    seg = size(popn,2) - (2 + ((VAR.transfers-1)*5));    % Disregarding the misc info, how many thrust/angle
                    seg = seg/VAR.transfers;                             % Thrust/angler per transfer
                    seg = seg/2;                                            % How many segments
                    
                    for i3 = 1:OPT_algo.Npop
                        
                        % From Departure to flyby 1
                        start = popn(i3,1);
                        tof = popn(i3,seg*2+5);
                        popn(i3,seg*2+6) = start+tof;
                    
                        % From flyby 2 to n
                        for i4 = 2:(VAR.transfers-1)
                            start = popn(i3,((i4-1)*((2*seg)+5))+1);
                            tof = popn(i3,i4*((2*seg)+5));
                            popn(i3,((i4)*((2*seg)+5))+1) = start+tof;
                        end
                        
                    end
                    
                end
                
            case {'PSO'}
                
                % Initials
                Vmax = OPT_algo.vmax*(VAR.high-VAR.low);

                % Pre-Allocate
                popn = zeros(OPT_algo.Npop,Nvar);
                bee(1:OPT_algo.Npop) = struct('pos',zeros(1,Nvar), 'vel',zeros(1,Nvar), 'f',0,...
                    'f_p',999999999999999, 'p',zeros(1,Nvar), 'f_g',999999999999999, 'g',zeros(1,Nvar)); 

                % Randomize the Values
                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        bee(i1).pos(i2) = randomNum(VAR.low(i2),VAR.high(i2),'dec');
                        bee(i1).vel(i2) = randomNum(-Vmax(i2),Vmax(i2),'dec');
                        if VAR.bin(i2)
                            bee(i1).pos(i2) = round(bee(i1).pos(i2));
                        end
                    end
                end
                
                
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end


                % Break the popn into transfer segments
                if VAR.transfers > 1
                    
                    % How many segments is the transfers broken into
                    seg = size(popn,2) - (2 + ((VAR.transfers-1)*5));    % Disregarding the misc info, how many thrust/angle
                    seg = seg/VAR.transfers;                             % Thrust/angler per transfer
                    seg = seg/2;                                            % How many segments
                    
                    for i3 = 1:OPT_algo.Npop
                        
                        % From Departure to flyby 1
                        start = bee(i3).pos(1);
                        tof = bee(i3).pos(seg*2+5);
                        bee(i3).pos(seg*2+6) = start+tof;
                    
                        % From flyby 2 to n
                        for i4 = 2:(VAR.transfers-1)
                            start = bee(i3).pos(((i4-1)*((2*seg)+5))+1);
                            tof = bee(i3).pos(i4*((2*seg)+5));
                            bee(i3).pos(((i4)*((2*seg)+5))+1) = start+tof;
                        end
                        popn(i3,:) = bee(i3).pos;
                    end
                    
                end
                
            case {'MBH'}
                
                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.N1_Outer,Nvar); %population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.N1_Outer
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end
                
                % Break the popn into transfer segments
                if VAR.transfers > 1
                    
                    % How many segments is the transfers broken into
                    seg = size(popn,2) - (2 + ((VAR.transfers-1)*5));    % Disregarding the misc info, how many thrust/angle
                    seg = seg/VAR.transfers;                             % Thrust/angler per transfer
                    seg = seg/2;                                            % How many segments
                    
                    for i3 = 1:OPT_algo.N1_Outer
                        
                        % From Departure to flyby 1
                        start = popn(i3,1);
                        tof = popn(i3,seg*2+5);
                        popn(i3,seg*2+6) = start+tof;
                    
                        % From flyby 2 to n
                        for i4 = 2:(VAR.transfers-1)
                            start = popn(i3,((i4-1)*((2*seg)+5))+1);
                            tof = popn(i3,i4*((2*seg)+5));
                            popn(i3,((i4)*((2*seg)+5))+1) = start+tof;
                        end
                        
                    end
                    
                end
                
            otherwise
                
                errorPathDisplay();
                errorAlgorithm();
                return
                
        end
        
    case {'MGALT_IN_FBSM_2D'}
        
        switch algorithm
            
            case {'GA','DE'}
                
                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.Npop,Nvar); %population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end
                
                % Break the popn into transfer segments
                for i3 = 1:OPT_algo.Npop
                    for i4 = 1:(VAR.transfers-1)
                        start = popn(i3,i4*11-10);
                        tof = popn(i3,i4*11);
                        popn(i3,i4*11+1) = start+tof;
                    end
                end
                
            case {'PSO'}
                
                % Initials
                Vmax = OPT_algo.vmax*(VAR.high-VAR.low);

                % Pre-Allocate
                popn = zeros(OPT_algo.Npop,Nvar);
                bee(1:OPT_algo.Npop) = struct('pos',zeros(1,Nvar), 'vel',zeros(1,Nvar), 'f',0,...
                    'f_p',999999999999999, 'p',zeros(1,Nvar), 'f_g',999999999999999, 'g',zeros(1,Nvar)); 

                % Randomize the Values
                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.Npop
                    for i2 = 1:Nvar
                        bee(i1).pos(i2) = randomNum(VAR.low(i2),VAR.high(i2),'dec');
                        bee(i1).vel(i2) = randomNum(-Vmax(i2),Vmax(i2),'dec');
                        if VAR.bin(i2)
                            bee(i1).pos(i2) = round(bee(i1).pos(i2));
                        end
                    end
                end

                % Break the popn into transfer segments
                % ie, popn is 75x40 with 2 transfers, so the first stansfer 
                % segment is popn(:,1:20) and the second segment is popn(:,21:end)
                for i3 = 1:OPT_algo.Npop
                    for i4 = 1:(VAR.transfers-1)
                        start = bee(i3).pos(i4*11-10);
                        tof = bee(i3).pos(i4*11);
                        bee(i3).pos(i4*11+1) = start+tof;
                    end
                    popn(i3,:) = bee(i3).pos;
                end
                
            case {'MBH'}
                
                % Initials
                bee = 0;

                % Pre-Allocate with Random Normalized Values
                popn = rand(OPT_algo.N1_Outer,Nvar); % Population

                % Adjust variables to be within bounds, round binary variables, set
                % constant variables equal to constant
                for i1 = 1:OPT_algo.N1_Outer
                    for i2 = 1:Nvar
                        if ~VAR.bin(i2)      % Non-binary variables
                            popn(i1,i2) = (VAR.high(i2)-VAR.low(i2))*popn(i1,i2) + VAR.low(i2);  
                        elseif VAR.bin(i2)   % Binary variables
                            popn(i1,i2) = round(popn(i1,i2));
                        end
                    end
                end
                
                % Break the popn into transfer segments
                for i3 = 1:OPT_algo.N1_Outer
                    for i4 = 1:(VAR.transfers-1)
                        start = popn(i3,i4*11-10);
                        tof = popn(i3,i4*11);
                        popn(i3,i4*11+1) = start+tof;
                    end
                end
                
            otherwise
                
                errorPathDisplay();
                errorAlgorithm();
                return
                
        end
        
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end



end


