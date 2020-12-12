function [] = plotOrbits(BOD,CONST,OPT,VAR,plot_vars)
% FORM: [] = plotOrbits(BOD,CONST,OPT,VAR,plot_vars)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Function designed to plot MGALT transfers
% |
% |     -This funciton plots the departure planet, flyby planet(s), and 
% |     the target planet along with the spacecraft's trajectories between 
% |     all of them. The funciton also plots the thruster pointing angle 
% |     for all of the transfers
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



%% User Choices

% Size of Sun and Planets
planet_size = 30;
sun_size = 40;

% General line sizes and shapes
line_width_orbit = 2;
line_width_thruster_angle = 2;
line_width_terminal_point = 5;
line_shape_terminal_point = 'kd';

% Direct method points and arrows
point_size = 25;
arrow_length = 0.3;
arrow_head_size = 1;
arrow_width = 2;

% Font
font_size = 24;
font_style = 'Times';
font_terminal = 'Terminal Points';

% Colors
sun_color = [0.95,0.95,0];
arrow_color = [0,0.5,0];



%% Constants

% Get Constants
mew_sun = CONST.Sun_mu;     % km^3/s^2
AU = CONST.AU;              % km/AU
SOI_divisions = 0:pi/180:2*pi;
orbit_margin = 1.01;

% Plot colors
plot_col = get(gca,'colororder');



%% Plot Orbits

switch OPT.solver
    
    case {'LT_DIR_FSM_2D'}

        %***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        % Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size);  % Plot the Sun
        hold on
        grid on
        Legend1{1} = 'Sun';
            
        % Plot the Departure Body Position
        body_data = plot_vars.planetary_conditions(:,1);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(1,:),'marker','.','markersize',planet_size);    % Body Position

        % Plot the Departure Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{1},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(1,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        % Plot the Departure Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{1},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,1),OPT.ode,mew_sun);
        body_orbit_line(1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');      % Exclude line from legend
        thrust_plot_dist = max([norm(body_orbit(1,:)),norm(body_orbit(2,:))]);                          	% Used for plotting, the length of the thrust pointing line
        clear body_period body_prop_period body_orbit

        % Plot the Spacecraft Transfer Arc
        sc_data = plot_vars.transfers{1,1};
        sc_pos_rad = sc_data(:,2);        	% DU
        sc_pos_ang = sc_data(:,3)*pi/180;	% rad
        [sc_X,sc_Y] = pol2cart(sc_pos_ang,sc_pos_rad);      % Convert to Cartesian
        plot(sc_X,sc_Y,'Color',plot_col(2,:),'linewidth',line_width_orbit)      % Spacecraft Trajectory

        % Calculate the Thrust Vectors
        thrust_switch = plot_vars.thrust_switch{1,1};                   % binary
        thrust_angle = plot_vars.thrust_phi{1,1};                       % deg
        thrust_R = plot_vars.seg_start{1,1}(:,1);                       % DU
        thrust_arc = plot_vars.seg_start{1,1}(:,2);                     % deg
        [thrust_X,thrust_Y] = pol2cart(thrust_arc*pi/180,thrust_R);     % points for the thrust vectors
        thrust_mag = arrow_length*ones(length(thrust_switch),1);       	% magnitude of thrust vectors for visual
        thrust_ang_global = thrust_arc + 90 - thrust_angle';         	% global angle for thrust vectors
        [thrust_X_global,thrust_Y_global] = pol2cart(thrust_ang_global*pi/180,thrust_mag);	% points for the thrust vectors 

        % Plot the Thrust Vectors
        quiver_count = 1;
        for i2 = 1:length(thrust_X)

            if thrust_switch(i2)
                thrust_point_line(1,i2) = plot(thrust_X(i2),thrust_Y(i2),...
                    'Color',arrow_color,'marker','.','markersize',point_size);
                thrust_quiver_line(1,quiver_count) = quiver(thrust_X(i2),thrust_Y(i2),thrust_X_global(i2),thrust_Y_global(i2),...
                    'color',[0,0.5,0],'MaxHeadSize',(arrow_head_size*AU),'LineWidth',arrow_width);
                set(get(get(thrust_quiver_line(1,quiver_count),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');	% Exclude line from legend
                quiver_count = quiver_count+1;
            else
                thrust_point_line(1,i2) = plot(thrust_X(i2),thrust_Y(i2),...
                    'Color',arrow_color,'marker','o','markersize',point_size/4);
            end

            set(get(get(thrust_point_line(1,i2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');              % Exclude line from legend

        end

        % Labels
        Legend1{2} = BOD.bodies{1};
        transfer_string1 = strcat("Transfer ",BOD.bodies{1}," to ",BOD.bodies{end});
        Legend1{3} = (transfer_string1);
        
        % Plot the Target Body
        body_data = plot_vars.planetary_conditions(:,end);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(2,:),'marker','.','markersize',planet_size);    % Body Position
        Legend1{4} = BOD.bodies{end};
        
        % Plot the Target Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{end},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        % Plot the Target Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{2},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),OPT.ode,mew_sun);
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
        
        % Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Direct FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        %***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        % Grid Basics
        hold on
        grid on
            
        % Plot Thruster Angle
        for i2 = 1:size(plot_vars.tspan{1,1},1)

            % If full or dashed
            if plot_vars.thrust_switch{1,1}(i2)
                thrust_segment_line(1,i2) = plot(plot_vars.tspan{1,1}(i2,:),...
                    [plot_vars.thrust_phi{1,1}(i2),plot_vars.thrust_phi{1,1}(i2)],...
                    'Color',plot_col(2,:),'linewidth',line_width_thruster_angle);
                set(get(get(thrust_segment_line(1,i2),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            else
                thrust_segment_line(1,i2) = plot(plot_vars.tspan{1,1}(i2,:),...
                    [plot_vars.thrust_phi{1,1}(i2),plot_vars.thrust_phi{1,1}(i2)],...
                    '--','Color',plot_col(2,:),'linewidth',line_width_thruster_angle);
                set(get(get(thrust_segment_line(1,i2),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            end

        end

        % Plot connecting lines
        for i3 = 2:size(plot_vars.tspan{1,1},1)

            % If full or dashed
            if plot_vars.thrust_switch{1,1}(i3)
                thrust_vertical_line(1,i3) = plot([plot_vars.tspan{1,1}(i3,1),plot_vars.tspan{1,1}(i3,1)],...
                    [plot_vars.thrust_phi{1,1}(i3-1),plot_vars.thrust_phi{1,1}(i3)],...
                    'Color',plot_col(2,:),'linewidth',line_width_thruster_angle);
                set(get(get(thrust_vertical_line(1,i3),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            else
                thrust_vertical_line(1,i3) = plot([plot_vars.tspan{1,1}(i3,1),plot_vars.tspan{1,1}(i3,1)],...
                    [plot_vars.thrust_phi{1,1}(i3-1),plot_vars.thrust_phi{1,1}(i3)],...
                    '--','Color',plot_col(2,:),'linewidth',line_width_thruster_angle);
                set(get(get(thrust_vertical_line(1,i3),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            end

        end

        % Plot Thruster Angle
        plot([0,0.00001],[0,00001],'Color',plot_col(2,:),'linewidth',line_width_thruster_angle);
        transfer_string2 = strcat("Transfer ",BOD.bodies{1}," to ",BOD.bodies{end});
        Legend2{1} = (transfer_string2);

        % Final Graph Elements
        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Direct FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
        
    case {'LT_IN_FSM_2D'}

        %***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        % Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size);  % Plot the Sun
        hold on
        grid on
        Legend1{1} = 'Sun';
            
        % Plot the Departure Body Position
        body_data = plot_vars.planetary_conditions(:,1);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(1,:),'marker','.','markersize',planet_size);    % Body Position

        % Plot the Departure Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{1},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(1,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        % Plot the Departure Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{1},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,1),OPT.ode,mew_sun);
        body_orbit_line(1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
        clear body_period body_prop_period body_orbit

        % Plot the Spacecraft Transfer Arc
        sc_data = plot_vars.transfers{1,1};
        sc_pos_rad = sc_data(:,5);            % DU
        sc_pos_ang = sc_data(:,6)*pi/180;   % rad;
        [sc_X,sc_Y] = pol2cart(sc_pos_ang,sc_pos_rad);    % Convert to Cartesian
        plot(sc_X,sc_Y,'Color',plot_col(2,:),'linewidth',line_width_orbit)	% Spacecraft Trajectory

        % Labels
        Legend1{2} = BOD.bodies{1};
        transfer_string1 = strcat("Transfer ",BOD.bodies{1}," to ",BOD.bodies{end});
        Legend1{3} = (transfer_string1);
        
        % Plot the Target Body
        body_data = plot_vars.planetary_conditions(:,end);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(2,:),'marker','.','markersize',planet_size);    % Body Position
        Legend1{4} = BOD.bodies{end};
        
        % Plot the Target Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{end},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        % Plot the Target Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{2},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),OPT.ode,mew_sun);
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
        
        % Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Indirect FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        %***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        % Grid Basics
        hold on
        grid on
        time = 0;
        phi = 0;

        % Thruster Data
        thruster_data = plot_vars.transfers{1,1};
        time = thruster_data(:,1)+time(end);
        lambda_1 = thruster_data(:,2);
        lambda_2 = thruster_data(:,3);

        % Run through all thruster angles
        for i2 = 1:length(lambda_1)
            phi(i2) = atan2d(-lambda_1(i2)/sqrt(lambda_1(i2)^2 + lambda_2(i2)^2),...
                            -lambda_2(i2)/sqrt(lambda_1(i2)^2 + lambda_2(i2)^2));
            if phi(i2) < 0
                phi(i2) = 360 + phi(i2);
            end
        end

        % Plot Thruster Angle
        plot(time,phi,'Color',plot_col(2,:),'linewidth',line_width_thruster_angle)
        transfer_string2 = strcat("Transfer ",BOD.bodies{1}," to ",BOD.bodies{end});
        Legend2{1} = (transfer_string2);
        
        % Final Graph Elements
        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Indirect FSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
    
    case {'MGALT_DIR_FBSM_2D'}
        
        %***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        % Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size);  % Plot the Sun
        hold on
        grid on
        count = 1;
        Legend1{count} = 'Sun';
        
        % Run through all the transfers
        for i1 = 1:VAR.transfers
            
            count = count+1;
            
            % Plot the Departure Body Position
            body_data = plot_vars.planetary_conditions(:,i1);
            body_X = body_data(1)/AU;	% DU
            body_Y = body_data(2)/AU; 	% DU
            plot(body_X,body_Y,'Color',plot_col(i1,:),'marker','.','markersize',planet_size);   % Body Position
            
            % Plot the Departure Body's SOI
            body_SOI = CONST.(strcat(BOD.bodies{i1},"_SOI"))/AU;
            body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
            body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
            body_SOI_line(i1,1) = plot(body_SOI_X,body_SOI_Y,'k--');
            set(get(get(body_SOI_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            % Plot the Departure Body's Orbit
            body_period = CONST.(strcat(BOD.bodies{i1},"_per"))*orbit_margin;
            body_prop_period = linspace(0,body_period*86400,500);
            [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,i1),OPT.ode,mew_sun);
            body_orbit_line(i1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
            set(get(get(body_orbit_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     % Exclude line from legend
            thrust_plot_dist = max([norm(body_orbit(1,:)),norm(body_orbit(2,:))]);                              % Used for plotting, the length of the thrust pointing line
            clear body_period body_prop_period body_orbit
            
            % Plot the Spacecraft Forward Transfer Arc
            sc_data_fs = plot_vars.transfers_fs{i1,1};
            sc_pos_rad_fs = sc_data_fs(:,2);            % DU
            sc_pos_ang_fs = sc_data_fs(:,3)*pi/180;     % rad
            [sc_X_fs,sc_Y_fs] = pol2cart(sc_pos_ang_fs,sc_pos_rad_fs);    % Convert to Cartesian
            plot(sc_X_fs,sc_Y_fs,'Color',plot_col(i1+1,:),'linewidth',line_width_orbit)	% Spacecraft Trajectory
            sc_fs_terminal_point(i1,1) = plot(sc_X_fs(end),sc_Y_fs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend

            % Plot the Spacecraft Backward Transfer Arc
            sc_data_bs = plot_vars.transfers_bs{i1,1};
            sc_pos_rad_bs = sc_data_bs(:,2);            % DU
            sc_pos_ang_bs = sc_data_bs(:,3)*pi/180;     % rad
            [sc_X_bs,sc_Y_bs] = pol2cart(sc_pos_ang_bs,sc_pos_rad_bs);    % Convert to Cartesian
            sc_transfer_line(i1,1) = plot(sc_X_bs,sc_Y_bs,'Color',plot_col(i1+1,:),'linewidth',line_width_orbit);	% Spacecraft Trajectory
            sc_bs_terminal_point(i1,1) = plot(sc_X_bs(end),sc_Y_bs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_transfer_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend


            % Calculate the Thrust Vectors
            thrust_switch_fs = plot_vars.thrust_switch_fs{i1,1};    % binary
            thrust_switch_bs = plot_vars.thrust_switch_bs{i1,1};
            thrust_angle_fs = plot_vars.thrust_phi_fs{i1,1};        % deg
            thrust_angle_bs = plot_vars.thrust_phi_bs{i1,1};
            thrust_R_fs = plot_vars.seg_start_fs{i1,1}(:,1);       	% DU
            thrust_R_bs = plot_vars.seg_start_bs{i1,1}(:,1);
            thrust_arc_fs = plot_vars.seg_start_fs{i1,1}(:,2);      % deg
            thrust_arc_bs = plot_vars.seg_start_bs{i1,1}(:,2);
            [thrust_X_fs,thrust_Y_fs] = pol2cart(thrust_arc_fs*pi/180,thrust_R_fs);	% points for the thrust vectors
            [thrust_X_bs,thrust_Y_bs] = pol2cart(thrust_arc_bs*pi/180,thrust_R_bs);
            thrust_mag_fs = arrow_length*ones(length(thrust_switch_fs),1);
            thrust_mag_bs = arrow_length*ones(length(thrust_switch_bs),1);
            thrust_ang_global_fs = thrust_arc_fs + 90 - thrust_angle_fs';        	% global angle for thrust vectors
            thrust_ang_global_bs = thrust_arc_bs + 90 - thrust_angle_bs';
            [thrust_X_global_fs,thrust_Y_global_fs] = pol2cart(thrust_ang_global_fs*pi/180,thrust_mag_fs);	% points for the thrust vectors 
            [thrust_X_global_bs,thrust_Y_global_bs] = pol2cart(thrust_ang_global_bs*pi/180,thrust_mag_bs);
            
            
            % Plot the Forward Shooting Thrust Vectors
            quiver_count_fs = 1;
            for i2 = 1:length(thrust_X_fs)
                
                if thrust_switch_fs(i2)
                    thrust_point_line_fs(i1,i2) = plot(thrust_X_fs(i2),thrust_Y_fs(i2),...
                        'Color',arrow_color,'marker','.','markersize',point_size);
                    thrust_quiver_line_fs(i1,quiver_count_fs) = ...
                        quiver(thrust_X_fs(i2),thrust_Y_fs(i2),...
                        thrust_X_global_fs(i2),thrust_Y_global_fs(i2),...
                        'color',[0,0.5,0],'MaxHeadSize',(arrow_head_size*AU),...
                        'LineWidth',arrow_width);
                    set(get(get(thrust_quiver_line_fs(i1,quiver_count_fs),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                              % Exclude line from legend
                    quiver_count_fs = quiver_count_fs+1;
                else
                    thrust_point_line_fs(i1,i2) = plot(thrust_X_fs(i2),thrust_Y_fs(i2),...
                        'Color',arrow_color,'marker','o','markersize',point_size/4);
                end
                
                set(get(get(thrust_point_line_fs(i1,i2),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');             % Exclude line from legend
            end
            
            % Plot the Backward Shooting Thrust Vectors
            quiver_count_bs = 1;
            for i3 = 1:length(thrust_X_bs)
                
                if thrust_switch_bs(i3)
                    thrust_point_line_bs(i1,i3) = plot(thrust_X_bs(i3),thrust_Y_bs(i3),...
                        'Color',arrow_color,'marker','.','markersize',point_size);
                    thrust_quiver_line_bs(i1,quiver_count_bs) = ...
                        quiver(thrust_X_bs(i3),thrust_Y_bs(i3),...
                        thrust_X_global_bs(i3),thrust_Y_global_bs(i3),...
                        'color',[0,0.5,0],'MaxHeadSize',(arrow_head_size*AU),...
                        'LineWidth',arrow_width);
                    set(get(get(thrust_quiver_line_bs(i1,quiver_count_bs),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');                              % Exclude line from legend
                    quiver_count_bs = quiver_count_bs+1;
                else
                    thrust_point_line_bs(i1,i3) = plot(thrust_X_bs(i3),thrust_Y_bs(i3),...
                        'Color',arrow_color,'marker','o','markersize',point_size/4);
                end
                
                set(get(get(thrust_point_line_bs(i1,i3),'Annotation'),...
                    'LegendInformation'),'IconDisplayStyle','off');             % Exclude line from legend
            end
            
            
            % Labels
            Legend1{count} = BOD.bodies{i1};
            transfer_string1 = strcat("Transfer ",BOD.bodies{i1}," to ",BOD.bodies{i1+1});
            count = count+1;
            Legend1{count} = (transfer_string1);
            
        end
        
        % Plot the Target Body
        body_data = plot_vars.planetary_conditions(:,end);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(i1+1,:),'marker','.','markersize',planet_size);    % Body Position
        Legend1{count+1} = BOD.bodies{end};
        
        % Plot the Target Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{end},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        % Plot the Target Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{end},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),OPT.ode,mew_sun);
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
        
        % Have 1 match point displayed
        sc_bs_terminal_point(i1,1) = plot(sc_X_bs(end),sc_Y_bs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
        Legend1{count+2} = font_terminal;
        
        % Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Direct FBSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        %***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        % Grid Basics
        hold on
        grid on
        time_fs = 0;
        time_bs = 0;
        count = 0;
        
        % Run through all the transfers
        for i1 = 1:VAR.transfers
            
            count = count+1;
            
            % Thruster Data
            time_fs = plot_vars.tspan_fs{i1,1} + time_bs(end,end);
            time_bs = plot_vars.tspan_bs{i1,1};
            time_bs = rot90(time_bs,2) + time_fs(end,end);
            
            
            % Plot Forward Thruster Angle
            for i2 = 1:size(plot_vars.tspan_fs{i1,1},1)
                
                % If full or dashed
                if plot_vars.thrust_switch_fs{i1,1}(i2)
                    thrust_segment_line_fs(i1,i2) = ...
                        plot([time_fs(i2,1),time_fs(i2,end)],...
                        [plot_vars.thrust_phi_fs{i1,1}(i2),...
                        plot_vars.thrust_phi_fs{i1,1}(i2)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_segment_line_fs(i1,i2),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                else
                    thrust_segment_line_fs(i1,i2) = ...
                        plot([time_fs(i2,1),time_fs(i2,end)],...
                        [plot_vars.thrust_phi_fs{i1,1}(i2),...
                        plot_vars.thrust_phi_fs{i1,1}(i2)],...
                        '--','Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_segment_line_fs(i1,i2),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                end
                
            end
            
            % Plot Forward Connecting Lines
            for i3 = 2:size(plot_vars.tspan_fs{i1,1},1)
                
                % If full or dashed
                if plot_vars.thrust_switch_fs{i1,1}(i3)
                    thrust_vertical_line_fs(i1,i3) = ...
                        plot([time_fs(i3,1),time_fs(i3,1)],...
                        [plot_vars.thrust_phi_fs{i1,1}(i3-1),...
                        plot_vars.thrust_phi_fs{i1,1}(i3)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_vertical_line_fs(i1,i3),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                else
                    thrust_vertical_line_fs(i1,i3) = ...
                        plot([time_fs(i3,1),time_fs(i3,1)],...
                        [plot_vars.thrust_phi_fs{i1,1}(i3-1),...
                        plot_vars.thrust_phi_fs{i1,1}(i3)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_vertical_line_fs(i1,i3),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                end
                
            end
            
            % Plot Backward Thruster Angle
            for i4 = 1:size(plot_vars.tspan_bs{i1,1},1)
                
                % If full or dashed
                if plot_vars.thrust_switch_bs{i1,1}(i4)
                    thrust_segment_line_bs(i1,i4) = ...
                        plot([time_bs(i4,1),time_bs(i4,end)],...
                        [plot_vars.thrust_phi_bs{i1,1}(i4),...
                        plot_vars.thrust_phi_bs{i1,1}(i4)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_segment_line_bs(i1,i4),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                else
                    thrust_segment_line_bs(i1,i4) = ...
                        plot([time_bs(i4,1),time_bs(i4,end)],...
                        [plot_vars.thrust_phi_bs{i1,1}(i4),...
                        plot_vars.thrust_phi_bs{i1,1}(i4)],...
                        '--','Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_segment_line_bs(i1,i4),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                end
                
            end
            
            % Plot Backward Connecting Lines
            for i5 = 2:size(plot_vars.tspan_bs{i1,1},1)
                
                % If full or dashed
                if plot_vars.thrust_switch_bs{i1,1}(i5)
                    thrust_vertical_line_bs(i1,i5) = ...
                        plot([time_bs(i5,1),time_bs(i5,1)],...
                        [plot_vars.thrust_phi_bs{i1,1}(i5-1),...
                        plot_vars.thrust_phi_bs{i1,1}(i5)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_vertical_line_bs(i1,i5),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                else
                    thrust_vertical_line_bs(i1,i5) = ...
                        plot([time_bs(i5,1),time_bs(i5,1)],...
                        [plot_vars.thrust_phi_bs{i1,1}(i5-1),...
                        plot_vars.thrust_phi_bs{i1,1}(i5)],...
                        'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
                    set(get(get(thrust_vertical_line_bs(i1,i5),...
                        'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off');     %Exclude line from legend
                end
                
            end
            
            
            % Plot Match Points
            sc_fs_terminal_point(i1,1) = plot(time_fs(i2,end),plot_vars.thrust_phi_fs{i1,1}(i2),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            sc_bs_terminal_point(i1,1) = plot(time_bs(1,1),plot_vars.thrust_phi_bs{i1,1}(1),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            
            
            % Plot Thruster Angle
            plot([0,0.00001],[0,00001],'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
            transfer_string2 = strcat("Transfer ",BOD.bodies{i1}," to ",BOD.bodies{i1+1});
            Legend2{count} = (transfer_string2);
            
        end

        % Final Graph Elements
        sc_bs_terminal_point(i1,1) = plot(time_bs(1,1),plot_vars.thrust_phi_bs{end,1}(1),line_shape_terminal_point,'linewidth',line_width_terminal_point);
        Legend2{count+1} = font_terminal;

        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Direct FBSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
    
    case {'MGALT_IN_FBSM_2D'}
    
        %***FIGURE 1: ORBITAL TRANSFERS***
        figure(1)
        
        % Grid Basics
        plot(0,0,'color',sun_color,'marker','.','markersize',sun_size);  % Plot the Sun
        hold on
        grid on
        count = 1;
        Legend1{count} = 'Sun';
        
        % Run through all the transfers
        for i1 = 1:VAR.transfers
            
            count = count+1;
            
            % Plot the Departure Body Position
            body_data = plot_vars.planetary_conditions(:,i1);
            body_X = body_data(1)/AU;	% DU
            body_Y = body_data(2)/AU; 	% DU
            plot(body_X,body_Y,'Color',plot_col(i1,:),'marker','.','markersize',planet_size);   % Body Position
            
            % Plot the Departure Body's SOI
            body_SOI = CONST.(strcat(BOD.bodies{i1},"_SOI"))/AU;
            body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
            body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
            body_SOI_line(i1,1) = plot(body_SOI_X,body_SOI_Y,'k--');
            set(get(get(body_SOI_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            % Plot the Departure Body's Orbit
            body_period = CONST.(strcat(BOD.bodies{i1},"_per"))*orbit_margin;
            body_prop_period = linspace(0,body_period*86400,500);
            [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,i1),OPT.ode,mew_sun);
            body_orbit_line(i1,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
            set(get(get(body_orbit_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            clear body_period body_prop_period body_orbit
            
            % Plot the Spacecraft Forward Transfer Arc
            sc_data_fs = plot_vars.transfers_fs{i1,1};
            sc_pos_rad_fs = sc_data_fs(:,5);            % DU
            sc_pos_ang_fs = sc_data_fs(:,6)*pi/180;   % rad;
            [sc_X_fs,sc_Y_fs] = pol2cart(sc_pos_ang_fs,sc_pos_rad_fs);    % Convert to Cartesian
            plot(sc_X_fs,sc_Y_fs,'Color',plot_col(i1+1,:),'linewidth',line_width_orbit)	% Spacecraft Trajectory
            sc_fs_terminal_point(i1,1) = plot(sc_X_fs(end),sc_Y_fs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            
            % Plot the Spacecraft Backwards Transfer Arc
            sc_data_bs = plot_vars.transfers_bs{i1,1};
            sc_pos_rad_bs = sc_data_bs(:,5);            % DU
            sc_pos_ang_bs = sc_data_bs(:,6)*pi/180;   % rad;
            [sc_X_bs,sc_Y_bs] = pol2cart(sc_pos_ang_bs,sc_pos_rad_bs);    % Convert to Cartesian
            sc_transfer_line(i1,1) = plot(sc_X_bs,sc_Y_bs,'Color',plot_col(i1+1,:),'linewidth',line_width_orbit);	% Spacecraft Trajectory
            sc_bs_terminal_point(i1,1) = plot(sc_X_bs(end),sc_Y_bs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_transfer_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            
            % Labels
            Legend1{count} = BOD.bodies{i1};
            transfer_string1 = strcat("Transfer ",BOD.bodies{i1}," to ",BOD.bodies{i1+1});
            count = count+1;
            Legend1{count} = (transfer_string1);
            
        end
        
        % Plot the Target Body
        body_data = plot_vars.planetary_conditions(:,end);
        body_X = body_data(1)/AU;	% DU
        body_Y = body_data(2)/AU; 	% DU
        plot(body_X,body_Y,'Color',plot_col(i1+1,:),'marker','.','markersize',planet_size);    % Body Position
        Legend1{count+1} = BOD.bodies{end};
        
        % Plot the Target Body's SOI
        body_SOI = CONST.(strcat(BOD.bodies{end},"_SOI"))/AU;
        body_SOI_X = (body_SOI/2)*cos(SOI_divisions) + body_X;
        body_SOI_Y = (body_SOI/2)*sin(SOI_divisions) + body_Y;
        body_SOI_line(end,1) = plot(body_SOI_X,body_SOI_Y,'k--');
        set(get(get(body_SOI_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        % Plot the Target Body's Orbit
        body_period = CONST.(strcat(BOD.bodies{end},"_per"))*orbit_margin;
        body_prop_period = linspace(0,body_period*86400,500);
        [~,body_orbit] = ode45(@orbit3D,body_prop_period,plot_vars.planetary_conditions(:,end),OPT.ode,mew_sun);
        body_orbit_line(end,1) = plot(body_orbit(:,1)/AU,body_orbit(:,2)/AU,'k--');
        set(get(get(body_orbit_line(end,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
        
        % Have 1 match point displayed
        sc_bs_terminal_point(i1,1) = plot(sc_X_bs(end),sc_Y_bs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
        Legend1{count+2} = font_terminal;
        
        % Final Graph Elements
        axis square
        hold off
        xlabel('X Position (AU)')
        ylabel('Y Positon (AU)')
        legend(string(Legend1),'location','eastoutside')
        title('Low-Thrust Orbit Transfer with Indirect FBSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        
        
        
        %***FIGURE 2: THRUSTER ANGLE***
        figure(2)
        
        % Grid Basics
        hold on
        grid on
        time_fs = 0;
        time_bs = 0;
        count = 0;
        
        % Run through all the transfers
        for i1 = 1:VAR.transfers
            
            count = count+1;
            
            % Thruster Data
            thruster_fs = plot_vars.transfers_fs{i1,1};
            time_fs = thruster_fs(:,1)+time_bs(end);
            lambda_1_fs = thruster_fs(:,2);
            lambda_2_fs = thruster_fs(:,3);
            
            thruster_bs = plot_vars.transfers_bs{i1,1};
            thruster_bs = flipud(thruster_bs);
            time_bs = thruster_bs(:,1)+time_fs(end);
            lambda_1_bs = thruster_bs(:,2);
            lambda_2_bs = thruster_bs(:,3);
            
            % Preallocate
            phi_fs = zeros(length(time_fs),1);
            phi_bs = zeros(length(time_bs),1);
            
            % Run through all thruster angles
            for i2 = 1:length(lambda_1_fs)
                phi_fs(i2) = atan2d(-lambda_1_fs(i2)/sqrt(lambda_1_fs(i2)^2 + lambda_2_fs(i2)^2),...
                    -lambda_2_fs(i2)/sqrt(lambda_1_fs(i2)^2 + lambda_2_fs(i2)^2));
                if phi_fs(i2) < 0
                    phi_fs(i2) = 360 + phi_fs(i2);
                end
            end
            
            for i3 = 1:length(lambda_1_bs)
                phi_bs(i3) = atan2d(-lambda_1_bs(i3)/sqrt(lambda_1_bs(i3)^2 + lambda_2_bs(i3)^2),...
                    -lambda_2_bs(i3)/sqrt(lambda_1_bs(i3)^2 + lambda_2_bs(i3)^2));
                if phi_bs(i3) < 0
                    phi_bs(i3) = 360 + phi_bs(i3);
                end
            end
            
            % Plot Thruster Angle
            plot(time_fs,phi_fs,'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
            sc_bs_line(i1,1) = plot(time_bs,phi_bs,'Color',plot_col(i1+1,:),'linewidth',line_width_thruster_angle);
            set(get(get(sc_bs_line(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            
            sc_fs_terminal_point(i1,1) = plot(time_fs(end),phi_fs(end),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_fs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            sc_bs_terminal_point(i1,1) = plot(time_bs(1),phi_bs(1),line_shape_terminal_point,'linewidth',line_width_terminal_point);
            set(get(get(sc_bs_terminal_point(i1,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');     %Exclude line from legend
            
            transfer_string2 = strcat("Transfer ",BOD.bodies{i1}," to ",BOD.bodies{i1+1});
            Legend2{count} = (transfer_string2);
            
        end
        
        % Final Graph Elements
        sc_bs_terminal_point(i1,1) = plot(time_bs(1),phi_bs(1),line_shape_terminal_point,'linewidth',line_width_terminal_point);
        Legend2{count+1} = font_terminal;
        
        hold off
        xlabel('Time (TU)')
        ylabel('\phi (deg)')
        legend(string(Legend2),'location','eastoutside')
        title('Thruster Pointing Angle \phi for Indirect FBSM')
        set(gca,'FontSize',font_size)
        set(gca,'FontName',font_style)
        ylim([0 360])
        yticks([0,45,90,135,180,225,270,315,360])
        
    otherwise
        
        errorPathDisplay();
        errorSolver();
        return
        
end



end


