function [isBroken] = errorMain(BOD,OPT)
% FORM: [isBroken] = errorMain(BOD,OPT)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function is called from the main script and will abort
% |     execution if any common user errors are input at the start
% |
% |     -The errors are currently:
% |         Variable limits
% |         MBH, percent elements ~= migration elements
% |         Old call to DIR_FSMseg_2D and forcing heliocentric orbits
% |         Optimzation before 1 Jan 1970 or after 31 Dec 2030
% |         
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
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -is_broken          (1,1)       [boolean]   	[unitless]
% |         If any of the user selections will cause the program to break
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initial

isBroken = false;



%% Variable Limits 

if size(BOD.bodies,2)-1 ~= size(OPT.thrust.tt_end,1)
    
    errorPathDisplay();
    fprintf(2,'The number of transfers and the list of bodies is not equal.\n')
    fprintf(2,'Fix "tof_total" or "tof_margin".\n')
    fprintf(2,'Fix the list of bodies under "bodies".\n')
    
    isBroken = true;
    return
    
end



%% MBH number of elements in percent selection ~= number of migrations

if isfield(OPT,'MBH')
    
    for i1 = 1:size(OPT.MBH,2)
    
        if (size(OPT.MBH(i1).per_feas,2) < (OPT.island.Nmig+1))...
                || (size(OPT.MBH(i1).per_rand,2) < (OPT.island.Nmig+1))

            errorPathDisplay();
            fprintf(2,'MBH will break on the last migration.\n')
            fprintf(2,'The number of migrations is larger than the size of "per_feas" or "per_rand" in "OPT.MBH(%1.0f)".\n',i1)
            fprintf(2,'Reduce the number of migrations or add another element to "per_feas" or "per_rand".\n')

            isBroken = true;
            return

        elseif (size(OPT.MBH(i1).per_feas,2) > (OPT.island.Nmig+1))...
                || (size(OPT.MBH(i1).per_rand,2) > (OPT.island.Nmig+1))

            errorPathDisplay();
            fprintf(2,'MBH has more elements in "per_feas" or "per_rand" than Island Iterations.\n')
            fprintf(2,'MBH has %2.0f "per_feas" parameters and %2.0f "per_rand" parameters.\n',size(OPT.MBH(i1).per_feas,2),size(OPT.MBH(i1).per_rand,2))
            fprintf(2,'The optimizer will still run but only the first %2.0f MBH parameter(s) will be used.\n\n',OPT.island.Nmig+1)
            fprintf(2,'Pausing for 10 sec...\n')
            pause(10)

        end
    
    end

end



%% Check OPT.thrust

% Having the segmented trajectory have heliocentric orbits within it
try
    
    if strcmp(OPT.thrust.orbit_check,'on')

        errorPathDisplay();
        fprintf(2,'As of right now, MGALT STOpS does not support FORCED multiple heliocentric revolutions.\n')
        fprintf(2,'This is future work, as the direct solver needs to specify which transfer arc to apply the heliocentric revolution to.\n')
        
        isBroken = true;
        return

    end

catch
    
    % For the indirect methods, which don't have the h-orbit field, make
    % sure to let the warning bypass the error catch
    warning('Reference to non-existent field')
    
end


% If doing MGALT Transfers
switch OPT.solver
    
    case {'MGALT_DIR_FBSM_2D','MGALT_IN_FBSM_2D'}
        
        % Check duty cycle
        if ~isfield(OPT.thrust,'duty_cycle')
            
            errorPathDisplay();
            fprintf(2,'For MGALT trajectories, the thruster duty cycle needs to be known.\n')
            
            isBroken = true;
            return
        end
        
        % Check number of thrusters
        if ~isfield(OPT.thrust,'n_available')
            
            errorPathDisplay();
            fprintf(2,'For MGALT trajectories, the total number of thrusters needs to be known.\n')
            
            isBroken = true;
            return
        end
        
    otherwise
        
        % Act as a pass
        
end


% Prevent the INDIRECT Method from having variable thrust
switch OPT.solver
    
    case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
        
        if strcmp(OPT.thrust.thrust_method,'variable')
            
            errorPathDisplay();
            fprintf(2,'The INDIRECT Method does not support varaible thrust.\n')
            
            isBroken = true;
            return
            
        end
        
    otherwise
        
        % Pass
    
end



%% Calculate a trajectory before 1 Jan 1970 or past 31 Dec 2030

JD1 = getJulianDate(BOD.window1); 
JD2 = getJulianDate(BOD.window2) + sum(OPT.thrust.tt_end) + ...
    sum(OPT.thrust.time(:,end));

if JD1 < 2440587.5
    
    errorPathDisplay();
    fprintf(2,'MGALT STOpS does not support trajectory optimization before 1 Jan 1970.\n')
    fprintf(2,'To negate this, new JPL Horizons Data is needed or accurate ode45 support for planet position will have to be implimented.\n')
    
    isBroken = true;
    return
    
end

if JD2 > 2462867.25
    
    errorPathDisplay();
    fprintf(2,'MGALT STOpS does not support trajectory optimization past 31 Dec 2030.\n')
    fprintf(2,'To negate this, new JPL Horizons Data is needed or accurate ode45 support for planet position will have to be implimented.\n')
    
    isBroken = true;
    return
    
end



%% Make sure Window2/JD2 ~< Window1/JD1

% JD's
if JD2 < JD1
    
    errorPathDisplay();
    fprintf(2,'Julian Date 2 before Julian Date 1, optimization will not work.\n')

    isBroken = true;
    return
                        
end

% Launch Windows
date1 = datenum(BOD.window1);
date2 = datenum(BOD.window2);
if date2 < date1
    
    errorPathDisplay();
    fprintf(2,'Launch window 2 is before launch window 1, optimization will not work.\n')

    isBroken = true;
    return
    
end



%% User Specified Algorithms ~= Algorithm List

% Get the island list and the number of those islands from Isl_opt
[GC,GR] = groupcounts(OPT.island.isl_list);

% Run through all of the island lists
for i1 = 1:size(GR,1)
    
    % For GA
    if strcmp(GR(i1,1),{'GA'})
        if GC(i1) < size(OPT.GA,2)
            fprintf('There are more GA islands from "parametersGA" than there are GA connections from "optionsIsland".\n')
            fprintf('The optimizer will still run but only the first %2.0f GA parameter(s) will be used.\n\n',GC(i1))
            fprintf(2,'Pausing for 10 sec...\n')
            pause(10)
        end
        if GC(i1) > size(OPT.GA,2)
            errorPathDisplay()
            fprintf(2,'There are more GA connections from "optionsIsland" than there are GA islands from "parametersGA".\n')
            isBroken = true;
            return
        end
    end
    
    % For DE
    if strcmp(GR(i1,1),{'DE'})
        if GC(i1) < size(OPT.DE,2)
            fprintf('There are more DE islands from "parametersDE" than there are DE connections from "optionsIsland".\n')
            fprintf('The optimizer will still run but only the first %2.0f DE parameter(s) will be used.\n\n',GC(i1))
            fprintf(2,'Pausing for 10 sec...\n')
            pause(10)
        end
        if GC(i1) > size(OPT.DE,2)
            errorPathDisplay()
            fprintf(2,'There are more DE connections from "optionsIsland" than there are DE islands from "parametersDE".\n')
            isBroken = true;
            return
        end
    end
    
    % For PSO
    if strcmp(GR(i1,1),{'PSO'})
        if GC(i1) < size(OPT.PSO,2)
            fprintf('There are more PSO islands from "parametersPSO" than there are PSO connections from "optionsIsland".\n')
            fprintf('The optimizer will still run but only the first %2.0f PSO parameter(s) will be used.\n\n',GC(i1))
            fprintf(2,'Pausing for 10 sec...\n')
            pause(10)
        end
        if GC(i1) > size(OPT.PSO,2)
            errorPathDisplay()
            fprintf(2,'There are more PSO connections from "optionsIsland" than there are PSO islands from "parametersPSO".\n')
            isBroken = true;
            return
        end
    end
    
    % For MBH
    if strcmp(GR(i1,1),{'MBH'})
        if GC(i1) < size(OPT.MBH,2)
            fprintf('There are more MBH islands from "parametersMBH" than there are MBH connections from "optionsIsland".\n')
            fprintf('The optimizer will still run but only the first %2.0f MBH parameter(s) will be used.\n\n',GC(i1))
            fprintf(2,'Pausing for 10 sec...\n')
            pause(10)
        end
        if GC(i1) > size(OPT.MBH,2)
            errorPathDisplay()
            fprintf(2,'There are more MBH connections from "optionsIsland" than there are MBH islands from "parametersMBH".\n')
            isBroken = true;
            return
        end
    end
       
end



%% Prevent Multiple FSM

switch OPT.solver
    
    case {'LT_DIR_FSM_2D','LT_IN_FSM_2D'}
        
        if size(BOD.bodies,2) > 2
            
            errorPathDisplay();
            fprintf(2,'MGALT STOpS does not support FSM with more than two planets.\n')
            fprintf(2,'It is recommended to use the FBSM instead.\n')

            isBroken = true;
            return
            
        end
        
    otherwise
   
        % Acts as a pass
        
end



%% Check if Parallel Toolbox

try 
    
	ver('distcomp');
    clc
    
catch
    
    errorPathDisplay();
    fprintf(2,'Parallel Processing Toolbox not installed.\n')
    fprintf(2,'Set "OPT.parallel" to false under the "Main Script" header.\n')
    fprintf(2,'If this message is an error, disable it under "errorMain".\n')

    isBroken = true;
    return
    
end



end


