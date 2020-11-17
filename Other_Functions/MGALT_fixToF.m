function [member] = MGALT_fixToF(BOD,OPT,VAR,member)
% FORM: [member] = MGALT_fixToF(BOD,OPT,VAR,member)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Function designed to fix the ToF constraints on the popn 
% |     variable. This method was originally done in the genPopn function, 
% |     but it was discovered that by running through the algorithms, 
% |     T1_JD+T1_ToF ~= T2_JD, leading to incirrect results where the 
% |     s/c would arrive at a transfers planet and then depart a few 
% |     weeks later. This funciton fixes this by running through the popn 
% |     and correcting any ToF errors
% |
% |     -This function needs to be called AFTER the popn is perterbated 
% |     in an algorithm funciton or else it won;t accomplish its purpose
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
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -member             (1,Nvar)    [float]         [unitless]
% |         A single member of a population input into the function
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -member             (1,Nvar)    [float]         [unitless]
% |         A single member of a population input into the function
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -Called from:
% |         MGALT_GA_mating
% |         MGALT_DE_nextGeneration
% |         MGALT_PSO_nextGeneration
% |         MGALT_MBH_randomize
% |
% |-----------------------------------------------------------------------



%% Fix the ToF

% For the different solver methods
switch OPT.solver
    
    case {'LT_IN_FSM_2D','LT_DIR_FSM_2D'}

        % Correction for transfer before ToF
        if member(1) < BOD.bodies_JD(1)
            diff_first = BOD.bodies_JD(1)-member(1);
            member(1) = member(1)+diff_first;
        end
        
        % Correction for the ToF being longer than acceptable
        if member(end) < (OPT.tof_total-OPT.tof_margin(1))
            member(end) = OPT.tof_total-OPT.tof_margin(1);
        end
        
        % Correction for the ToF being longer than acceptable
        if member(end) > (OPT.tof_total+OPT.tof_margin(2))
            member(end) = OPT.tof_total+OPT.tof_margin(2);
        end
        
    case {'MGALT_IN_FBSM_2D'}
        
        % Correction for first transfer before ToF
        if member(1) < BOD.bodies_JD(1)
            diff_first = BOD.bodies_JD(1)-member(1);
            member(1) = member(1)+diff_first;
        end
        
        switch VAR.transfers
            
            case {1}

                % Correction for the ToF being past the JD
                if (member(1)+member(8))> BOD.bodies_JD(end)
                   diff_end = BOD.bodies_JD(end) - (member(1)+member(8));
                   member(8) = diff_end;
                end
                
            otherwise
                
                % Corrections for ToF between
                for i1 = 1:(VAR.transfers)-1
                    start = member(i1*11-10);
                    tof = member(i1*11);
                    member(i1*11+1) = start+tof;
                end
                
                % Correction for the ToF being past the JD
                last_JD = member(end-7);
                last_ToF = member(end);
                
                if (last_JD+last_ToF) > BOD.bodies_JD(end)
                    diff_end = BOD.bodies_JD(end) - (last_JD+last_ToF);
                    for i2 = 1:VAR.transfers
                        member(i2*11-10) = member(i2*11-10) + diff_end;
                    end
                    
                end
                
        end
        
    case {'MGALT_DIR_FBSM_2D'}
        
        % Correction for first transfer before ToF
        if member(1) < BOD.bodies_JD(1)
            diff_first = BOD.bodies_JD(1)-member(1);
            member(1) = member(1)+diff_first;
        end
        
        switch VAR.transfers
            
            case {1}

                % Correction for the ToF being longer than acceptable
                if member(end) < (OPT.tof_total-OPT.tof_margin(1))
                    member(end) = OPT.tof_total-OPT.tof_margin(1);
                end

                % Correction for the ToF being longer than acceptable
                if member(end) > (OPT.tof_total+OPT.tof_margin(2))
                    member(end) = OPT.tof_total+OPT.tof_margin(2);
                end
                
            otherwise
                
                % How many segments is the transfers broken into
                seg = size(member,2) - (2 + ((VAR.transfers-1)*5));    % Disregarding the misc info, how many thrust/angle
                seg = seg/VAR.transfers;                             % Thrust/angler per transfer
                seg = seg/2; 
                
                % Corrections for ToF between
                start = member(1);
                tof = member(seg*2+5);
                member(seg*2+6) = start+tof;
                for i1 = 2:VAR.transfers-1
                    start = member(((i1-1)*((2*seg)+5))+1);
                    tof = member(i1*((2*seg)+5));
                    member(((i1)*((2*seg)+5))+1) = start+tof;
                end
                
                % Correction for the ToF being past the JD
                last_JD = member((end-1)-(2*seg));
                last_ToF = member(end);
                
                if (last_JD+last_ToF) > BOD.bodies_JD(end)
                    diff_end = BOD.bodies_JD(end) - (last_JD+last_ToF);
                    
                    member(1) = member(1)+diff_end;
                    for i2 = 2:VAR.transfers-1
                        member(((i2-1)*((2*seg)+5))+1) = member(((i2-1)*((2*seg)+5))+1)+diff_end;
                    end
                    member((end-1)-(2*seg)) = member((end-1)-(2*seg))+diff_end;
                    
                end
                
        end
        
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end



end


