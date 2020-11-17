function [R,V] = MGALT_stateBodies(JD_target,BOD,CONST,OPT,VAR,array_bodies)
% FORM: [R,V] = MGALT_stateBodies(JD_target,BOD,CONST,OPT,VAR,array_bodies)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function to extract planetary location and velocity for certain 
% |     Julan Dates from large variable arrays
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -JD_target        	(1,1)       [float]         [JD]
% |         The desired JD to get planetary R and V at
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
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -array_bodies       (3,1)       [int]        [unitless]
% |         An array containing the index numbers for the bodies
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -R                  (6,1)       [float]         [AU]
% |         Radius vector [planet1; planet2]
% |     -V                  (6,1)       [float]         [AU/TU]
% |         Velocity vector [planet1; planet2]
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Get Info

% Find which indicies the 0 value is located at
[~,index] = min(abs(BOD.bodies_JD-JD_target));
is_backwards = (BOD.bodies_JD(index)-JD_target) >= 0;
JD_delta = BOD.bodies_JD(index)-JD_target;          % How many JDs ahead of the target JD the closest segment is 

% Check if JD_delta*86400 > -1 and < 1, due to tspan being indexed by 1
if ((JD_delta*86400 > -1) && (JD_delta*86400 < 1))
    
    R = BOD.bodies_R(array_bodies,index);
    V = BOD.bodies_V(array_bodies,index);
    return
    
end

% Check error for first or last index under/overflow
if ((index == 1) && (JD_delta >= 1)) || ((index == size(BOD.bodies_JD,2)) && (JD_delta <= 1))
    
    errorPathDisplay();
    fprintf(2,'Error with JD parsing!\n')
    fprintf(2,'The desired JD date is out of range of the valid JD dates stored in "BOD.bodies_JD".\n')
    fprintf(2,'ODE45 will make solution inaccurate...exiting program.\n\n')
    
    return
        
end

% If is_backwards, means backwards propogate to get to R and V values
if is_backwards
    tspan = 0:-1:(-JD_delta*86400);
else
    tspan = 0:(-JD_delta*86400);
end

% Get R and V data
R_init = BOD.bodies_R(array_bodies,index);
V_init = BOD.bodies_V(array_bodies,index);

% Use ODE45 to get the planet position at the desired JD
[~,data] = ode45(@orbit3D,tspan,[R_init;V_init],OPT.ode,CONST.Sun_mu);

% Final Results
R = data(end,1:3)';
V = data(end,4:6)';



end


