function [] = displayResults(BOD,CONST,OPT,VAR,run_time,eval_info)
% FORM: [] = displayResults(BOD,CONST,OPT,VAR,run_time,eval_info)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function designed to display optimization results
% |
% |     -This funciton gives info regarding optimization time and 
% |     iterations, orbital transfers, and other information regarding 
% |     some of the inputs used for STOpS
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
% |     -run_time        	(1,1)       [float]        	[min]
% |         How long the optimization process took to run
% |     -eval_info          (1,1)       [struct]        [unitless]
% |         Best solutions, best costs, number of iterations, etc...
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Plot Results

%Section Description:
%{
-Check through all of the results in order to fund the best results for the
  optimization methods.
%}


% ***Calculated***
final_soln = zeros(OPT.island.Nisl,size(VAR.bin,2));
final_f = zeros(1,OPT.island.Nisl);

for i3 = 1:OPT.island.Nisl
    final_soln(i3,:) = eval_info(OPT.island.Nmig+1,i3).optimal_soln(1,:);
    final_f(i3) = eval_info(OPT.island.Nmig+1,i3).f_best(i3);
end

[~,I] = sort(final_f);
sorted_final_soln = final_soln(I,:);

% Plot the results
[f,plot_vars] = feval(OPT.solver,sorted_final_soln(1,:),BOD,CONST,OPT,VAR);
plotOrbits(BOD,CONST,OPT,VAR,plot_vars)



%% Display Completion Message

clc
fprintf('[\b---------------------]\b\n')
fprintf('[\bOptimization Complete]\b\n')
fprintf('[\b---------------------]\b\n\n')
fprintf('   Solver Used: %s\n\n',OPT.solver)
fprintf('Total Run Time: %5.2f minutes\n',run_time)
fprintf('                %3.2f hours\n',run_time/60)
fprintf('      Parallel: %1.0f \n',OPT.parallel)

total = 0;
for i1 = 1:size(eval_info,1)
    for i2 = 1:size(eval_info,2)
        total = total+eval_info(i1,i2).total_evals;
    end
end

fprintf('   Total Iters: %5.0f\n',total);
fprintf('     Best Cost: %7.4f\n',f)
fprintf('\n\n\n')



%% Orbit Information

TU = CONST.TU;

fprintf('[\b-----------------]\b\n')
fprintf('[\bOrbit Information]\b\n')
fprintf('[\b-----------------]\b\n\n')
fprintf('                *  DD/MM/YYYY  HH:MM:SS  *\n\n')

[~,~,~,~,~,~,~,date] = julian2greg(sorted_final_soln(1,1));
fprintf('   Departure Body: %s\n',BOD.bodies{1,1})
fprintf('   Departure Date: %2.0f/%2.0f/%2.0f, %2.0f:%2.0f:%2.0f\n',date(1),date(2),date(3),date(4),date(5),date(6))
fprintf('\n====>\n\n')


% Display the transfers
switch VAR.transfers
    
    case {1}
        
        % Time fo Flight
        switch OPT.solver
            case {'LT_DIR_FSM_2D'}
                ToF = plot_vars.tspan{1,1}(end,end)*TU;
            case {'MGALT_DIR_FBSM_2D'}
                ToF = plot_vars.tspan_fs{1,1}(end,end)*2*TU;
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                ToF = plot_vars.tspan(end,end)*TU;
        end

        % Arrival Date
        [~,~,~,~,~,~,~,date] = julian2greg(plot_vars.JD(end,end));

        % Fuel Usage
        data = plot_vars.transfers{end,1};
        fuel_start = data(1,end);
        fuel_end = data(end,end);

        fprintf('     Arrival Body: %s\n',BOD.bodies{1,end})
        fprintf('   Time of Flight: %5.2f days\n',ToF)
        fprintf('     Arrival Date: %2.0f/%2.0f/%2.0f, %2.0f:%2.0f:%2.0f\n',date(1),date(2),date(3),date(4),date(5),date(6))
        fprintf('  Fuel Mass Start: %5.2f kg\n',fuel_start)
        fprintf('    Fuel Mass End: %5.2f kg\n',fuel_end)
        fprintf('\n***\n\n')

        % Fuel warning
        if fuel_end <= 1
            fprintf(2,'\nWARNING!\n')
            fprintf(2,'Fuel dropped below 1kg on Transfer Segment %1.0f\n',1)
            fprintf(2,'Most likely the orbit did not converge with %s\n',BOD.bodies{1,end})
            fprintf(2,'Re-run the solver with different engine parameters!\n')
        end

        fprintf('\n. . . . . . . . . . . . . . . .\n\n')
        
    otherwise
        
        for i4 = 2:VAR.transfers

            % Time fo Flight
            switch OPT.solver
                case {'LT_DIR_FSM_2D'}
                    ToF(i4-1) = plot_vars.tspan{i4-1,end}(end,end)*TU;
                case {'MGALT_DIR_FBSM_2D'}
                    ToF(i4-1) = plot_vars.tspan_fs{i4-1,1}(end,end)*2*TU;
                case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                    ToF(i4-1) = plot_vars.tspan(i4-1,end)*TU;
            end

            % Flyby Dates
            flyby = plot_vars.JD(i4-1,end);
            [~,~,~,~,~,~,~,date] = julian2greg(flyby);

            % Fuel Usage
            data = plot_vars.transfers{i4-1,1};
            fuel_start = data(1,end);
            fuel_end = data(end,end);

            fprintf('  Transfer Body %1.0f: %s\n',i4-1,BOD.bodies{1,i4})
            fprintf('   Time of Flight: %5.2f days\n',ToF(i4-1))
            fprintf('       Flyby Date: %2.0f/%2.0f/%2.0f, %2.0f:%2.0f:%2.0f\n\n',date(1),date(2),date(3),date(4),date(5),date(6))
            fprintf('  Fuel Mass Start: %5.2f kg\n',fuel_start)
            fprintf('    Fuel Mass End: %5.2f kg\n',fuel_end)

            % Fuel warning
            if fuel_end <= 0.5
                fprintf(2,'\nWARNING!\n')
                fprintf(2,'Fuel dropped below 0.5kg on Transfer Segment %1.0f\n',i4-1)
                fprintf(2,'Most likely the orbit did not converge with %s\n',BOD.bodies{1,end})
                fprintf(2,'Re-run the solver with different engine parameters!\n\n')
                fprintf(2,'This is just a warning and not an actual error, nothing was suspended during operation.\n')
            end

            fprintf('\n====>\n\n')

        end

        % Time fo Flight
        switch OPT.solver
            case {'LT_DIR_FSM_2D'}
                ToF(end+1) = plot_vars.tspan{end,end}(end,end)*TU;
            case {'MGALT_DIR_FBSM_2D'}
                ToF(end+1) = plot_vars.tspan_fs{end,end}(end,end)*2*TU;
            case {'LT_IN_FSM_2D','MGALT_IN_FBSM_2D'}
                ToF(end+1) = plot_vars.tspan(end,end)*TU;
        end

        % Arrival Date
        [~,~,~,~,~,~,~,date] = julian2greg(plot_vars.JD(end,end));

        % Fuel Usage
        data = plot_vars.transfers{end,1};
        fuel_start = data(1,end);
        fuel_end = data(end,end);

        fprintf('     Arrival Body: %s\n',BOD.bodies{1,end})
        fprintf('   Time of Flight: %5.2f days\n',ToF(end))
        fprintf('     Arrival Date: %2.0f/%2.0f/%2.0f, %2.0f:%2.0f:%2.0f\n\n',date(1),date(2),date(3),date(4),date(5),date(6))
        fprintf('  Fuel Mass Start: %5.2f kg\n',fuel_start)
        fprintf('    Fuel Mass End: %5.2f kg\n',fuel_end)

        % Fuel warning
        if fuel_end <= 1
            fprintf(2,'\nWARNING!\n')
            fprintf(2,'Fuel dropped below 1kg on Transfer Segment %1.0f\n',i4)
            fprintf(2,'Most likely the orbit did not converge with %s\n',BOD.bodies{1,end})
            fprintf(2,'Re-run the solver with different engine parameters!\n')
        end

        fprintf('\n. . . . . . . . . . . . . . . .\n\n')
        
end

% Print off the total flight time
fprintf('Total Flight Time: %5.2f days\n',sum(ToF))
fprintf('                 : %2.2f years\n\n',sum(ToF)/365.25)

%Print off the fuel mass fraction
mass_total = plot_vars.transfers{1,1}(1,end);
mass_structure = plot_vars.transfers{end,1}(end,end);
mass_fuel = mass_total-mass_structure;
mass_fraction = mass_fuel/mass_total * 100;
fprintf('   Prop mass frac: %2.2f percent\n',mass_fraction)

fprintf('\n\n\n')



%% Other Information

fprintf('[\b-----------------]\b\n')
fprintf('[\bOther Information]\b\n')
fprintf('[\b-----------------]\b\n\n')

% Search Window
fprintf('    Search Window:\n')
fprintf('    %4.0f %4.0f %4.0f\n',BOD.window1(1),BOD.window1(2),BOD.window1(3))
fprintf('    %4.0f %4.0f %4.0f\n',BOD.window2(1),BOD.window2(2),BOD.window2(3))
fprintf('\n. . . . . . . . . . . . . . . .\n\n')

% ToF Parameters
fprintf('    Desired ToF:\n\n')
for i5 = 1:VAR.transfers
   fprintf('    %s to %s\n',BOD.bodies{1,i5},BOD.bodies{1,i5+1})
   fprintf('    %4.0f Days\n',OPT.thrust.tt_end(i5))
   fprintf('    -%4.0f +%4.0f Days\n',OPT.thrust.time(i5,1),OPT.thrust.time(i5,2));
   fprintf('\n')
end
fprintf('\n. . . . . . . . . . . . . . . .\n\n')

% Thrust Parameters
fprintf('    Engine Parameters:\n')
fprintf('    %s\n',OPT.thrust.thrust_method)
fprintf('    @ %2.4f N\n',OPT.thrust.thrust)
fprintf('    @ %4.2f Isp\n',OPT.thrust.Isp)
fprintf('    @ %2.4e kg/s\n',OPT.thrust.mdot)
fprintf('\n. . . . . . . . . . . . . . . .\n\n')



end


