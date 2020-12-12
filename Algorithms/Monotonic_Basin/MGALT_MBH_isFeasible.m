function [feasible,pert,transfer] = ...
    MGALT_MBH_isFeasible(BOD,CONST,OPT,OPT_algo,VAR,plot_vars,per_feas)
% FORM: [feasible,pert,transfer] = ...
%       MGALT_MBH_isFeasible(BOD,CONST,OPT,VAR,plot_vars,per_feas)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function looks at all of the transfers from the previous 
% |     iteration of the members within "MGALT_MBH_function". If the 
% |     transfers do not intersect the planet, then they are determined 
% |     to not be feasible. If the transfers do intersect the planet, 
% |     then the whole member is considered to be feasible.
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
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -plot_vars          (1,1)       [struct]     	[unitless]
% |         An object containing a lot of information about the 
% |         optimization parameters including: transfers(t and y ode 
% |         outputs), thrust values, thruster pointing angles, transfer 
% |         starting position, planet start/end locations for each 
% |         transfer, JD of each transfer, and tspans of each transfer
% |     -per_feas           (1,1)       [float]     	[unitless]
% |         The current feasibility percentage
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -feasible           (1,1)       [boolean]   	[unitless]
% |         The new perturbed member consisting of the original member and
% |         the perturbations added to it
% |     -pert  	(inputs.transfers,4)   	[float]       	[JD]
% |         The pert locations for this MBH struct
% |     -transfer        	(1,1)       [int]           [unitless]
% |         The current transfer number
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -For future improvement, this function should look at individual 
% |     segments of the transfer(s). If an individual segment is feasible, 
% |     then it should get flagged as a potential candidate and cycled 
% |     through into the cluster. The segments of the member which are 
% |     not feasible should then be altered until they are feasible.
% |
% |     -Some of the above infrastructure is setup with the pert and
% |     transfer returns, but they were not fully implimented.
% |
% |-----------------------------------------------------------------------



%% Initials

feasible = false;
transfer = zeros(1,VAR.transfers);
count = 1;

% Get some constants
AU = CONST.AU;           %km/AU



%% Disregard Departure and ToF Members

% Compare against the times to see if random perturbations made the
% transfer outside of the defined limit
leave = zeros(VAR.transfers,2);
tof = zeros(VAR.transfers,2);

switch OPT.solver
    
    case {'LT_IN_FSM_2D','LT_DIR_FSM_2D'}
        
        % Run through all of the transfers
        for i1 = 1:VAR.transfers

            % Departure bounds
            leave(1,1) = plot_vars.JD(1,1) < VAR.low(1);
            leave(1,2) = plot_vars.JD(1,1) > VAR.high(1);

            % TOF Bounds
            tof(1,1) = plot_vars.JD(1,2) < (VAR.low(1)...
                + VAR.low(end));
            tof(1,2) = plot_vars.JD(1,2) > (VAR.high(1)...
                + VAR.high(end));

        end
        
    case {'MGALT_IN_FBSM_2D','MGALT_DIR_FBSM_2D'}
        
        % Departure bounds
        leave(1,1) = plot_vars.JD(1,1) < VAR.low(1);
        leave(1,2) = plot_vars.JD(1,1) > VAR.high(1);
        
        switch VAR.transfers
            
            case {1}
                
                tof(1,1) = plot_vars.JD(1,2) < (VAR.low(1) + VAR.low(end));
                tof(1,2) = plot_vars.JD(1,2) > (VAR.high(1) + VAR.high(end));
                
            otherwise
                
                tof(1,1) = plot_vars.JD(1,2) < (VAR.low(1) + VAR.low(11));
                tof(1,2) = plot_vars.JD(1,2) > (VAR.high(1) + VAR.high(11));
                
                % Run through all of the transfers
                for i1 = 1:(VAR.transfers-2)
                    leave(i1+1,1) = plot_vars.JD(i1+1,1) < (VAR.low(i1*11+1));
                    leave(i1+1,2) = plot_vars.JD(i1+1,1) > (VAR.high(i1*11+1));
                    tof(i1+1,1) = plot_vars.JD(i1+1,2) < (VAR.low(i1*11+1) + VAR.low(i1*11+11));
                    tof(i1+1,2) = plot_vars.JD(i1+1,2) > (VAR.high(i1*11+1) + VAR.high(i1*11+11));
                end
                
                leave(end,1) = plot_vars.JD(end,1) < VAR.low(end-7);
                leave(end,2) = plot_vars.JD(end,1) > VAR.high(end-7);
                tof(end,1) = plot_vars.JD(end,2) < (VAR.low(end-7) + VAR.low(end));
                tof(end,2) = plot_vars.JD(end,2) > (VAR.high(end-7) + VAR.high(end));
                
        end
        
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end

pert = [leave,tof];

% The perturbations will break the ToF, exit
if any(pert,'all')
    return
end



%% Compare Feasibilities

switch OPT.solver
    
    case {'LT_IN_FSM_2D','LT_DIR_FSM_2D'}	% Feasibility is defined as intersecting with the planet's SOI

        % Check to see if intersects with SOI
        [does_intersect,pos_body,values] = ...
        SOIIntersectCheck(BOD,CONST,OPT,VAR,plot_vars,2);

        % If it does intersect with the SOI
        if does_intersect
            SOI_rad = CONST.(strcat(BOD.bodies{2},"_SOI"))/AU/2;	% km/AU
            SOI_percent = ((SOI_rad*per_feas) + SOI_rad);         	% km/AU

            pos_sc = values(1,1:2);

            % Check to see if the terminal point is located within the
            % region defined as feasible
            if ( (norm(pos_sc) >= (norm(pos_body(1:2))-SOI_percent)) ...
                    && (norm(pos_sc) <= (norm(pos_body(1:2))+SOI_percent)) )

                feasible(1) = true;
                transfer(1) = 1;

            end

        end
        
    case {'MGALT_IN_FBSM_2D'}            	% Feasibility is defined as the patch points being within a certain tolerance of each other
        
        % Run through all of the transfers
        for i2 = 1:VAR.transfers
            
            forward = plot_vars.transfers_fs{i2,1};
            backward = plot_vars.transfers_bs{i2,1};
            
            % I need to see if the match points are feasible.
            % Dr. Englander recommended using 1e-5 AU for the tolerance
            % point
            tol = OPT_algo.feas_tol;
            check = zeros(1,3);
            
            % Get the position norms, vel rad, and vel tan
            norm_pos_fs = norm(forward(end,5),forward(end,6));
            norm_pos_bs = norm(backward(end,5),backward(end,6));
            
            vel_rad_fs = forward(end,7);
            vel_rad_bs = backward(end,7);
            
            vel_tan_fs = forward(end,8);
            vel_tan_bs = backward(end,8);
            
            % Check feasibility of position X and Y
            if (abs(norm_pos_fs-norm_pos_bs) <= (tol*per_feas))
                
                % Make sure in the same quadrant
                if (sign(forward(end,5)) == sign(backward(end,5))) && ...
                        (sign(forward(end,6)) == sign(backward(end,6)))
                    check(1) = 1;
                end
                
            end
            
            % Check feasibility of radial velocity
          	if (abs(vel_rad_fs-vel_rad_bs) <= (tol*per_feas))
                check(2) = 1;
            end
            
            % Check feasibility of tangential velocity
            if (abs(vel_tan_fs-vel_tan_bs) <= (tol*per_feas))
                check(3) = 1;
            end
            
            % If all of the checks pass
            if all(check)
                feasible(count) = true;
                transfer(count) = i2;
            else
                feasible(count) = false;
                transfer(count) = i2;
            end
            count = count+1;
             
        end
        
    case {'MGALT_DIR_FBSM_2D'}            	% Feasibility is defined as the patch points being within a certain tolerance of each other
        
        % Run through all of the transfers
        for i2 = 1:VAR.transfers
            
            forward = plot_vars.transfers_fs{i2,1};
            backward = plot_vars.transfers_bs{i2,1};
            
            % I need to see if the match points are feasible.
            % Dr. Englander recommended using 1e-5 AU for the tolerance
            % point
            tol = OPT_algo.feas_tol;
            check = zeros(1,3);
            
            % Get the position norms, vel rad, and vel tan
            norm_pos_fs = norm(forward(end,2),forward(end,3));
            norm_pos_bs = norm(backward(end,2),backward(end,3));
            
            vel_rad_fs = forward(end,4);
            vel_rad_bs = backward(end,4);
            
            vel_tan_fs = forward(end,5);
            vel_tan_bs = backward(end,5);
            
            % Check feasibility of position X and Y
            if (abs(norm_pos_fs-norm_pos_bs) <= (tol*per_feas))
                
                % Make sure in the same quadrant
                if (sign(forward(end,2)) == sign(backward(end,2))) && ...
                        (sign(forward(end,3)) == sign(backward(end,3)))
                    check(1) = 1;
                end
                
            end
            
            % Check feasibility of radial velocity
          	if (abs(vel_rad_fs-vel_rad_bs) <= (tol*per_feas))
                check(2) = 1;
            end
            
            % Check feasibility of tangential velocity
            if (abs(vel_tan_fs-vel_tan_bs) <= (tol*per_feas))
                check(3) = 1;
            end
            
            % If all of the checks pass
            if all(check)
                feasible(count) = true;
                transfer(count) = i2;
            else
                feasible(count) = false;
                transfer(count) = i2;
            end
            count = count+1;
            
        end
 
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end



end


