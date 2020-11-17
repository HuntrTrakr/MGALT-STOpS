function [times,thrust,phi] = DIR_sepVariables(solver,member,index,seg)
% FORM: [times,thrust,phi] = DIR_sepVariables(solver,member,index,seg)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function designed to parse time, thrust, and thrust pointing 
% |     angle arrays from the MGALT_DIR_FSM_2D method
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -solver             (1,n)       [string]        [unitless]
% |         The solver (cost function) handle which is used to determine 
% |         the correct solver to use. This string is used in an ode45 
% |     -member             (1,Nvar)    [float]         [unitless]
% |         A single member of a population input into the function
% |     -index              (1,1)       [int]       	[unitless]
% |         The index value of a for loop used as an input to other 
% |         functions. The index value is useful in etracting information 
% |         like the planet string, R/V arrays, etc... without having to 
% |         pass in extra/unuseful information
% |     -seg             	(1,1)       [int]           [unitless]
% |         The segments for the direct method 
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -times              (1,2)       [float]         [JD][Day]
% |         The JD and ToF parsed from the member array
% |     -thrust             (1,Nseg)	[bool][float] 	[unitless][N]
% |         Binary indicator if the s/c is thrusting or not
% |         The thrust for each angle
% |     -phi                (1,Nseg) 	[float]         [deg]
% |         Thrust poointing angle for each segment
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Sep Vars

switch solver
    
    case {'LT_DIR_FSM_2D'}

        times(1) = member(1);
        times(2) = member(end);

        count = 1;
        for i1 = 2:2:(size(member,2)-1)
            thrust(count) = member(i1);   % N
            phi(count) = member(i1+1);    % deg
            count = count + 1;
        end
        
    case {'MGALT_DIR_FBSM_2D'}
        
        % Only call this in the for i1 = 2:transfers-1 portion, as the
        % start and end locations are different
        times(1) = member(((index-1)*((2*seg)+5))+1);
        times(2) = member(index*((2*seg)+5));
        data = member(   (((index-1)*((2*seg)+5))+2) : (index*((2*seg)+5)-4)   );
        thrust = data(1:2:end-1);
        phi = data(2:2:end);
        
    otherwise
        
        errorPathDisplay()
        errorSolver()
        return
        
end
        


end


