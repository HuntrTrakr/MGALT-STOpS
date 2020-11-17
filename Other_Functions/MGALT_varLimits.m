function [info] = MGALT_varLimits(BOD,OPT)
% FORM: [info] = MGALT_varLimits(BOD,OPT)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function designed to generate variable limits and binary check 
% |     array for low-thrust trajectory variable strings
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
% |     -info               (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

% opt_thrust = OPT.thrust;
bodies = BOD.bodies;

% Desired transfer times for each segment
Htime = OPT.thrust.tt_end;
transfers = size(bodies,2)-1;



%% Method Selection

switch OPT.solver
    
    case {'LT_DIR_FSM_2D'}	% Segmented Method
        
        % Preallocate: Lows/Highs/Binary
        low = zeros(transfers,(OPT.thrust.Nseg*2)+2);
        high = zeros(transfers,(OPT.thrust.Nseg*2)+2);
        bin = zeros(transfers,(OPT.thrust.Nseg*2)+2);
        
        % Departure body
        low(1,1) = getJulianDate(BOD.window1);      % earliest departure time (JD)
        high(1,1) = getJulianDate(BOD.window2);     % latest departure time (JD)
        
        switch OPT.thrust.thrust_method
        
            case {'constant','equation'}
                
                % This thrust method assumes that "first Tswitch assumed 1"
                % from Sheehan's code base. It was kept the same way while
                % updating to MGALT, but only for the departure planet.
                % It was assumed only for the departure planet because
                % subsequent flybys may need the thrust to be off after
                % happening
                
                % Departure body
                low(1,2) = 1;   % Tswitch low (on)
                high(1,2) = 1;  % Tswitch high (on)
                bin(1,2) = 0;   % non-binary variable, to force the thrust to always be on for departure

                % Segments between departure and next planet
                for i1 = 3:2:OPT.thrust.Nseg*2+1     
                    low(1,i1) = 0;      % phi low
                    high(1,i1) = 360;   % phi high
                    bin(1,i1) = 0;      % non-binary variable
                    
                    low(1,i1+1) = 0;  	% Tswitch low (off)
                    high(1,i1+1) = 1; 	% Tswitch high (on)
                    bin(1,i1+1) = 1;   	% binary variable  
                end
                
                low(1,end) = Htime(1) - OPT.thrust.time(1,1);	% shortest transfer time (days)
                high(1,end) = Htime(1) + OPT.thrust.time(1,2);	% longest transfer time (days)
                bin(1,end) = 0;                                 % non-binary variable
                
                % Other transfer segment departures
                for i2 = 2:transfers
                    
                    low(i2,:) = low(1,:);
                    high(i2,:) = high(1,:);
                    bin(i2,:) = bin(1,:);
                    
                    low(i2,1) = low(i2-1,1)+low(i2-1,end);              % earliest departure time
                    high(i2,1) = high(i2-1,1)+high(i2-1,end);           % latest departure time
                    low(i2,end) = Htime(i2) - OPT.thrust.time(i2,1);	% shortest transfer time (days)
                    high(i2,end) = Htime(i2) + OPT.thrust.time(i2,2);	% longest transfer time (days)
                    
                end
                
                % Make the lows for every transfer be 0 so it's not the
                % assumption of having a powered thrust like from the
                % departure body
                low(2:end,2) = 0;
                bin(2:end,2) = 1;   % For having preliminary segments be powered/unpowered after transfer
                
            case {'variable'}
                
                % Departure body
                for i1 = 1:2:OPT.thrust.Nseg*2
                    low(1,i1+1) = OPT.thrust.thrust(1);     % T low
                    high(1,i1+1) = OPT.thrust.thrust(2);    % T high
                    
                    low(1,i1+2) = 0;                        % phi low
                    high(1,i1+2) = 360;                     % phi high
                end
                
                low(1,end) = Htime(1) - OPT.thrust.time(1,1);   % shortest transfer time (days)
                high(1,end) = Htime(1) + OPT.thrust.time(1,2);  % longest transfer time (days)
                bin(1,end) = 0;                                 % non-binary variable
                
                %Other transfer segment departures
                for i2 = 2:transfers
                    
                    low(i2,:) = low(1,:);
                    high(i2,:) = high(1,:);
                    
                    low(i2,1) = low(i2-1,1)+low(i2-1,end);      % earliest departure time
                    high(i2,1) = high(i2-1,1)+high(i2-1,end);	% latest departure time
                    low(i2,end) = Htime(i2) - OPT.thrust.time(i2,1);    % shortest transfer time (days)
                    high(i2,end) = Htime(i2) + OPT.thrust.time(i2,2);	% longest transfer time (days)
                    
                end
                
            otherwise
                
                errorPathDisplay();
                disp('ERROR')
                return
        
        end
        
    case {'LT_IN_FSM_2D'}    	% Costate Method
        
        % Preallocate: Lows/Highs/Binary
        low = zeros(transfers,5);
        high = zeros(transfers,5);
        bin = zeros(transfers,5);

        % Departure body
        low(1,1) = getJulianDate(BOD.window1);          % earliest departure time (JD)
        high(1,1) = getJulianDate(BOD.window2);         % latest departure time (JD)
        low(1,2:4) = -1;                                % minimum of -1 for lamda variables
        high(1,2:4) = 1;                                % maximium of 1 for lamda variables
        low(1,5) = Htime(1) - OPT.thrust.time(1,1);     % shortest transfer time (days)
        high(1,5) = Htime(1) + OPT.thrust.time(1,2);	% longest transfer time (days) 
        
        %Other transfer segment departures
        for i1 = 2:transfers
            low(i1,1) = low(i1-1,1)+low(i1-1,5);            % earliest departure time (JD)
            high(i1,1) = high(i1-1,1)+high(i1-1,5);         % latest departure time (JD)
            low(i1,2:4) = -1;                               % minimum of -1 for lamda variables
            high(i1,2:4) = 1;                               % maximium of 1 for lamda variables
            low(i1,5) = Htime(i1) - OPT.thrust.time(i1,1);  % shortest transfer time (days)
            high(i1,5) = Htime(i1) + OPT.thrust.time(i1,2); % longest transfer time (days) 
        end    
        
    case {'MGALT_DIR_FBSM_2D'}	% Segmented Method
        
        % Make sure segments are even
        if mod(OPT.thrust.Nseg,2)
            
            errorPathDisplay();
            fprintf(2,'The DIRECT FBSM only supports an even number of segments\n\n');
            return
        end
        
        % Preallocate: Lows/Highs/Binary
        low_dep = zeros(1,OPT.thrust.Nseg+1);
        high_dep = zeros(1,OPT.thrust.Nseg+1);
        bin_dep = zeros(1,OPT.thrust.Nseg+1);
        
        low_ari = zeros(1,OPT.thrust.Nseg+1);
        high_ari = zeros(1,OPT.thrust.Nseg+1);
        bin_ari = zeros(1,OPT.thrust.Nseg+1);

        low_trans = zeros(transfers-1,2*OPT.thrust.Nseg+5);
        high_trans = zeros(transfers-1,2*OPT.thrust.Nseg+5);
        bin_trans = zeros(transfers-1,2*OPT.thrust.Nseg+5);

        reshaped_low_trans = [];
        reshaped_high_trans = [];
        reshaped_bin_trans = [];

        % Departure body
        low_dep(1,1) = getJulianDate(BOD.window1);    	% earliest departure time (JD)
        high_dep(1,1) = getJulianDate(BOD.window2);     % latest departure time (JD)
        bin_dep(1,1) = 0;                               % non-binary variable
        
        
        % Thrust methods        
        switch OPT.thrust.thrust_method
        
            case {'constant','equation'}
                
                % This thrust method assumes that "first Tswitch assumed 1"
                % from Sheehan's code base. It was kept the same way while
                % updating to MGALT, but only for the departure planet.
                % It was assumed only for the departure planet because
                % subsequent flybys may need the thrust to be off after
                % happening
                
                % Departure body
                low_dep(1,2) = 1;   % Tswitch low (on)
                high_dep(1,2) = 1;  % Tswitch high (on)
                bin_dep(1,2) = 0;   % non-binary variable, to force the thrust to always be on for departure

                % Departure segments
                for i1 = 3:2:OPT.thrust.Nseg
                    low_dep(1,i1) = 0;      % phi low
                    high_dep(1,i1) = 360;   % phi high
                    
                    low_dep(1,i1+1) = 0;  	% Tswitch low (off)
                    high_dep(1,i1+1) = 1; 	% Tswitch high (on)
                    bin_dep(1,i1+1) = 1;   	% binary variable  
                end
                low_dep(1,OPT.thrust.Nseg+1) = 0;
                high_dep(1,OPT.thrust.Nseg+1) = 360;
                bin_dep(1,OPT.thrust.Nseg+1) = 0;
                
                % Arrival segments
                for i2 = 1:2:OPT.thrust.Nseg
                    low_ari(1,i2) = 0;      % Tswitch low (off)
                    high_ari(1,i2) = 1;     % Tswitch high (on)
                    bin_ari(1,i2) = 1;      % binary variable 
                    
                    low_ari(1,i2+1) = 0;    % phi low
                    high_ari(1,i2+1) = 360; % phi high
                end
                low_ari(1,OPT.thrust.Nseg+1) = Htime(end) - OPT.thrust.time(end,1);     % earliest departure time
                high_ari(1,OPT.thrust.Nseg+1) = Htime(end) + OPT.thrust.time(end,2);	% latest departure time
                bin_ari(1,OPT.thrust.Nseg+1) = 0;                                       % non-binary variable
                
                
                % For transfer segments
                if transfers > 1
                    
                    % All transfers
                    for i3 = 1:transfers-1
                        
                        % Prior to Flyby
                        for i4 = 1:2:OPT.thrust.Nseg
                            low_trans(i3,i4) = 0;       % Tswitch low (off)
                            high_trans(i3,i4) = 1;      % Tswitch high (on)
                            bin_trans(i3,i4) = 1;       % binary variable
                            
                            low_trans(i3,i4+1) = 0;    	% phi low
                            high_trans(i3,i4+1) = 360;  % phi high                
                        end
                        
                        % Flyby Variables
                        low_trans(i3,OPT.thrust.Nseg+1) = 0;                        % minimum of 0 for the flyby coefficient for rp
                        high_trans(i3,OPT.thrust.Nseg+1) = 1;                       % maximum of 1 for the flyby coefficient for rp
                        low_trans(i3,OPT.thrust.Nseg+2:OPT.thrust.Nseg+3) = -1;     % Englander up to unit vector control min
                        high_trans(i3,OPT.thrust.Nseg+2:OPT.thrust.Nseg+3) = 1;     % Englander up to unit vector control max
                        low_trans(i3,OPT.thrust.Nseg+4) = Htime(i3) - OPT.thrust.time(i3,1);                        % shortest transfer time (days)
                        high_trans(i3,OPT.thrust.Nseg+4) = Htime(i3) + OPT.thrust.time(i3,2);                       % longest transfer time (days) 
                        low_trans(i3,OPT.thrust.Nseg+5) = low_dep(1,1) + sum(low_trans(:,OPT.thrust.Nseg+4));       % earliest departure time (JD)
                        high_trans(i3,OPT.thrust.Nseg+5) = high_dep(1,1) + sum(high_trans(:,OPT.thrust.Nseg+4));	% latest departure time (JD)
                        
                        % After Flyby
                        for i5 = (OPT.thrust.Nseg+6):2:(2*OPT.thrust.Nseg+4)
                            low_trans(i3,i5) = 0;       % Tswitch low (off)
                            high_trans(i3,i5) = 1;      % Tswitch high (on)
                            bin_trans(i3,i5) = 1;       % binary variable
                            
                            low_trans(i3,i5+1) = 0;    	% phi low
                            high_trans(i3,i5+1) = 360;	% phi high
                        end

                    end
                    
                    % Reshape the vectors          
                    for i6 = 1:transfers-1
                        reshaped_low_trans = [reshaped_low_trans, low_trans(i6,:)];
                        reshaped_high_trans = [reshaped_high_trans, high_trans(i6,:)];
                        reshaped_bin_trans = [reshaped_bin_trans, bin_trans(i6,:)];
                    end
              
                end
                 
            case {'variable'}
                
                % This thrust method assumes that "first Tswitch assumed 1"
                % from Sheehan's code base. It was kept the same way while
                % updating to MGALT, but only for the departure planet.
                % It was assumed only for the departure planet because
                % subsequent flybys may need the thrust to be off after
                % happening
                
                % Departure body
                low_dep(1,2) = OPT.thrust.thrust(1);   	% T low
                high_dep(1,2) = OPT.thrust.thrust(2);   % T high
                bin_dep = zeros(1,OPT.thrust.Nseg+1);   % non-binary variable

                % Departure segments
                for i1 = 3:2:OPT.thrust.Nseg
                    
                    low_dep(1,i1) = 0;      % phi low
                    high_dep(1,i1) = 360;   % phi high
                    
                    low_dep(1,i1+1) = OPT.thrust.thrust(1);     % T low
                    high_dep(1,i1+1) = OPT.thrust.thrust(2);	% T high
                    
                end
                low_dep(1,OPT.thrust.Nseg+1) = 0;
                high_dep(1,OPT.thrust.Nseg+1) = 360;
                
                % Arrival segments
                for i2 = 1:2:OPT.thrust.Nseg
                    low_ari(1,i2) = OPT.thrust.thrust(1);	% T low
                    high_ari(1,i2) = OPT.thrust.thrust(2);	% T high
                    
                    low_ari(1,i2+1) = 0;    % phi low
                    high_ari(1,i2+1) = 360; % phi high
                end
                low_ari(1,OPT.thrust.Nseg+1) = Htime(end) - OPT.thrust.time(end,1);     % earliest departure time
                high_ari(1,OPT.thrust.Nseg+1) = Htime(end) + OPT.thrust.time(end,2);	% latest departure time
                
                
                % For transfer segments
                if transfers > 1
                    
                    % All transfers
                    for i3 = 1:transfers-1
                        
                        % Prior to Flyby
                        for i4 = 1:2:OPT.thrust.Nseg
                            low_trans(i3,i4) = OPT.thrust.thrust(1);	% T low
                            high_trans(i3,i4) = OPT.thrust.thrust(2);	% T high
                            
                            low_trans(i3,i4+1) = 0;    	% phi low
                            high_trans(i3,i4+1) = 360;  % phi high                
                        end
                        
                        % Flyby Variables
                        low_trans(i3,OPT.thrust.Nseg+1) = 0;                        % minimum of 0 for the flyby coefficient for rp
                        high_trans(i3,OPT.thrust.Nseg+1) = 1;                       % maximum of 1 for the flyby coefficient for rp
                        low_trans(i3,OPT.thrust.Nseg+2:OPT.thrust.Nseg+3) = -1;     % Englander up to unit vector control min
                        high_trans(i3,OPT.thrust.Nseg+2:OPT.thrust.Nseg+3) = 1;     % Englander up to unit vector control max
                        low_trans(i3,OPT.thrust.Nseg+4) = Htime(i3) - OPT.thrust.time(i3,1);                        % shortest transfer time (days)
                        high_trans(i3,OPT.thrust.Nseg+4) = Htime(i3) + OPT.thrust.time(i3,2);                       % longest transfer time (days) 
                        low_trans(i3,OPT.thrust.Nseg+5) = low_dep(1,1) + sum(low_trans(:,OPT.thrust.Nseg+4));       % earliest departure time (JD)
                        high_trans(i3,OPT.thrust.Nseg+5) = high_dep(1,1) + sum(high_trans(:,OPT.thrust.Nseg+4));	% latest departure time (JD)
                        
                        % After Flyby
                        for i5 = (OPT.thrust.Nseg+6):2:(2*OPT.thrust.Nseg+4)
                            low_trans(i3,i5) = OPT.thrust.thrust(1);	% T low
                            high_trans(i3,i5) = OPT.thrust.thrust(2);	% T high
                            
                            low_trans(i3,i5+1) = 0;    	% phi low
                            high_trans(i3,i5+1) = 360;	% phi high
                        end

                    end
                    
                    % Reshape the vectors          
                    for i6 = 1:transfers-1
                        reshaped_low_trans = [reshaped_low_trans, low_trans(i6,:)];
                        reshaped_high_trans = [reshaped_high_trans, high_trans(i6,:)];
                        reshaped_bin_trans = [reshaped_bin_trans, bin_trans(i6,:)];
                    end
              
                end
                
            otherwise
                
                errorPathDisplay();
                disp('ERROR')
                return
        
        end
        
	case {'MGALT_IN_FBSM_2D'}    	% Costate Method
        
        % Preallocate: Lows/Highs/Binary
        low_trans = zeros(transfers-1,11);
        high_trans = zeros(transfers-1,11);
        bin_trans = zeros(transfers-1,11);

        reshaped_low_trans = [];
        reshaped_high_trans = [];
        reshaped_bin_trans = [];

        
        % Departure body
        low_dep(1,1) = getJulianDate(BOD.window1);    	% earliest departure time (JD)
        high_dep(1,1) = getJulianDate(BOD.window2);     % latest departure time (JD)
        low_dep(1,2:4) = -1;                         	% minimum of -1 for dep lamda variables
        high_dep(1,2:4) = 1;                          	% maximium of 1 for dep lamda variables
        bin_dep(1,1:4) = 0;                             % non-binary variable
        
        % Arrival Body
        low_ari(1,1:3) = -1;                                    % minimum of -1 for ari lamda variables
        high_ari(1,1:3) = 1;                                    % maximium of 1 for ari lamda variables
        low_ari(1,4) = Htime(end) - OPT.thrust.time(end,1);     % shortest transfer time (days)
        high_ari(1,4) = Htime(end) + OPT.thrust.time(end,2);    % longest transfer time (days)
        bin_ari(1,1:4) = 0;                                     % non-binary variable
        
        
        % If there are transfers
        if transfers > 1
            
            % Transfer segments
            for i1 = 1:transfers-1
                bin_trans(i1,1:11) = 0;     % non-binary variable
                low_trans(i1,1:3) = -1;     % minimum of -1 for forward lamda variables
                high_trans(i1,1:3) = 1;     % maximum of 1 for forward lamda variables
                low_trans(i1,4) = 0;        % minimum of 0 for the flyby coefficient for rp
                high_trans(i1,4) = 1;       % maximum of 1 for the flyby coefficient for rp
                low_trans(i1,5:6) = -1;     % Englander up to unit vector control min
                high_trans(i1,5:6) = 1;     % Englander up to unit vector control max
                low_trans(i1,7) = Htime(i1) - OPT.thrust.time(i1,1);    % shortest transfer time (days)
                high_trans(i1,7) = Htime(i1) + OPT.thrust.time(i1,2);   % longest transfer time (days) 
                low_trans(i1,8) = low_dep(1,1)+sum(low_trans(:,7));         % earliest departure time (JD)
                high_trans(i1,8) = high_dep(1,1)+sum(high_trans(:,7));      % latest departure time (JD)
                low_trans(i1,9:11) = -1;     % minimum of -1 for backward lamda variables
                high_trans(i1,9:11) = 1;     % maximum of 1 for backward lamda variables
            end

            % Reshape the vectors          
            for i2 = 1:transfers-1
                reshaped_low_trans = [reshaped_low_trans, low_trans(i2,:)];
                reshaped_high_trans = [reshaped_high_trans, high_trans(i2,:)];
                reshaped_bin_trans = [reshaped_bin_trans, bin_trans(i2,:)];
            end
        
        end
    
    otherwise                               % Error
        
        errorPathDisplay();
        errorSolver();
        return
        
end



%% Reshape the size

switch OPT.solver
    
    case {'LT_DIR_FSM_2D','LT_IN_FSM_2D'}
        
        reshaped_low(1:size(bin,2)) = low(1,:);
        reshaped_high(1:size(bin,2)) = high(1,:);
        reshaped_bin(1:size(bin,2)) = bin(1,:);

        for i1 = 2:transfers
            reshaped_low = [reshaped_low, low(i1,:)];
            reshaped_high = [reshaped_high, high(i1,:)];
            reshaped_bin = [reshaped_bin, bin(i1,:)];
        end
        
    case {'MGALT_DIR_FBSM_2D','MGALT_IN_FBSM_2D'}
        
        reshaped_low = [low_dep, reshaped_low_trans, low_ari];
        reshaped_high = [high_dep, reshaped_high_trans, high_ari];
        reshaped_bin = [bin_dep, reshaped_bin_trans, bin_ari];
        
end

info.low = reshaped_low;
info.high = reshaped_high;
info.bin = reshaped_bin;
info.transfers = transfers;



end


