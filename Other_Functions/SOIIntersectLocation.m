function [does_intersect,return_values] = ...
    SOIIntersectLocation(BOD,CONST,body_pos,sc_pos,sc_vel,index)
% FORM: [does_intersect,return_values] = ...
%       SOIIntersectLocation(BOD,CONST,body_pos,sc_pos,sc_vel,index)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function performs two different operations depending on what 
% |     results are found within the provided data. The first operation 
% |     takes the input planetary location, the planet, and the whole s/c 
% |     transfer arc. It then checks every pair of points along the 
% |     transfer arc and sees if it intersects with the planet's SOI. 
% |     If no intersection is detected, the function will exit.
% |
% |     -If SOI intersection is detected, the function will determine the 
% |     heliocentric X and Y, as well as vX and vY locations of the SOI 
% |     intersect. This is done by identifying the intersect segment and 
% |     finding out exactly where the segment intersects the circle which 
% |     is the SOI. From there, a ratio of linear interpolation is taken 
% |     to extract the vX and vY values.
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
% |     -body_pos        	(1,1)       [float]         [AU]
% |         Heliocentric X,Y components of the desired transfer body at a
% |         predetermined JD
% |     -sc_pos             (n,3)       [float]         [AU]
% |         Heliocentric Radius components for the s/c transfer
% |     -sc_vel             (n,3)       [float]       	[AU/TU]
% |         Heliocentric Velocity components for the s/c transfer
% |     -index              (1,1)       [int]       	[unitless]
% |         Index for the current body number
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -does_intersect 	(1,1)       [boolean]    	[unitless]
% |         If the s/c trajectory intersects the target SOI at any point
% |         along the orbital path
% |     -return_values  	(2,2)       [float]         [AU][AU/TU]
% |         Position (X,Y); Velocity (X,Y) for the spacecraft
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -Math taken from
% |         http://www.ambrsoft.com/TrigoCalc/Circles2/circlrLine_.htm
% |         https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
% |
% |     -Future work, located at line 86
% |         Run through the entirity of the body location through the 
% |         orbit, not just the final point...this would require a lot of 
% |         restructuring with the ode terminitation if a planet position 
% |         was located
% |
% |-----------------------------------------------------------------------



%% Initials

% Constants
AU = CONST.AU;                                              % km/AU
SOI_rad = CONST.(strcat(BOD.bodies{index},"_SOI"))/AU/2;	% km/AU

% Default returns
does_intersect = false;
return_values = NaN;

% Body position in canonical units
body_X = body_pos(1)/AU;
body_Y = body_pos(2)/AU;
% body_Z = body_pos(3)/AU;

% S/C position in canonical units
sc_X = sc_pos(:,1);
sc_Y = sc_pos(:,2);
% sc_Z = sc_pos(:,3);

% S/C velocity in canonical units
sc_vX = sc_vel(:,1);
sc_vY = sc_vel(:,2);
% sc_vZ = sc_vel(:,3);



%% Check

% *** FUTURE WORK ***
% for i1 = 1:size(body_X,1)

    % Where the centerpoint is located, this is the planet position
    center = [body_X, body_Y];

    % Run through the entirity of the s/c transfer
    for i2 = 1:size(sc_X,1)-1
        
        % The two points which define a line segment
        p0 = [sc_X(i2), sc_Y(i2)];      % Point 1 of the line
        p1 = [sc_X(i2+1), sc_Y(i2+1)];  % Point 2 of the line
        
        % Take the two points and get into slope intercept form
        % y = mx+b
        line_m = (p1(2)-p0(2))/(p1(1)-p0(1));
        line_d = p0(2) - line_m*p0(1);
        
        % Discriminant
        % This will allow knowledge if the line segment will intersect
        % with the circle...if true then more math to calculate if it
        % actually happens over the segment or if in infinity
        discriminant = (SOI_rad^2)*(1 + line_m^2) - ((center(2) - line_m*center(1) - line_d)^2);
        
        % Line which will not intersect
        if discriminant < 0 || isnan(discriminant)
            continue;   % Goes back to the start of the for loop above
        end

        
        % Now move onto the line segment which will intersect
        % Will need to calculate if it actually passes through the SOI
        % radius
        length_segment = sqrt((p1(1)-p0(1))^2 + (p1(2)-p0(2))^2);
        
        % Intersection points X and Y
        inter_x1 = (center(1) + center(2)*line_m - line_d*line_m + sqrt(discriminant))/(1 + line_m^2);
        inter_x2 = (center(1) + center(2)*line_m - line_d*line_m - sqrt(discriminant))/(1 + line_m^2);
        inter_y1 = (line_d + center(1)*line_m + center(2)*line_m^2 + line_m*sqrt(discriminant))/(1+line_m^2);
        inter_y2 = (line_d + center(1)*line_m + center(2)*line_m^2 - line_m*sqrt(discriminant))/(1+line_m^2);

        % Intersection points in vector form
        inter_1_pos = [inter_x1,inter_y1];
        inter_2_pos = [inter_x2,inter_y2];
        
        % Make a vector of the length of the intersection points to compare
        % to, see if the intersection points lie on the line segment
        length_inter(1) = sqrt((p0(1)-inter_1_pos(1))^2 + (p0(2)-inter_1_pos(2))^2);
        length_inter(2) = sqrt((p0(1)-inter_2_pos(1))^2 + (p0(2)-inter_2_pos(2))^2);
        length_inter(3) = sqrt((p1(1)-inter_1_pos(1))^2 + (p1(2)-inter_1_pos(2))^2);
        length_inter(4) = sqrt((p1(1)-inter_2_pos(1))^2 + (p1(2)-inter_2_pos(2))^2);
        
        % Run through the vector and check if the intersection points are
        % on the line segment
        flag = zeros(1,size(length_inter,2));
        for i3 = 1:length(length_inter)
            if length_inter(i3) <= length_segment
                flag(i3) = true;
            else
                flag(i3) = false;
            end
        end
        
        % Check which intersection point is satisfied
        if not(flag)
            continue
            
        elseif (flag(1) && flag(3))     % Trajectory starts outside of SOI and enters into it
            does_intersect = true;
            
            % Get the linear interpolated velocity values
            [~,~,ratio] = linearInterpolation(p0,p1,inter_1_pos(1));
            pv0 = [sc_vX(i2), sc_vY(i2)];      % Point 1 of the line
            pv1 = [sc_vX(i2+1), sc_vY(i2+1)];  % Point 2 of the line
            [vX_new,vY_new,~] = linearInterpolation(pv0,pv1,NaN,ratio);
            inter_1_vel = [vX_new,vY_new];
            
            return_values = [inter_1_pos;inter_1_vel];
            return

        elseif (flag(2) && flag(4))     % Trajectory starts within the SOI and leaves it
            does_intersect = true;
            
            % Get the linear interpolated velocity values
            [~,~,ratio] = linearInterpolation(p0,p1,inter_2_pos(1));
            pv0 = [sc_vX(i2), sc_vY(i2)];      % Point 1 of the line
            pv1 = [sc_vX(i2+1), sc_vY(i2+1)];  % Point 2 of the line
            [vX_new,vY_new,~] = linearInterpolation(pv0,pv1,NaN,ratio);
            inter_2_vel = [vX_new,vY_new];
            
            return_values = [inter_2_pos;inter_2_vel];
            return
            
        elseif flag                     % Trajectory enters and exits the SOI...find the direction of the transfer
            does_intersect = true;
            
            
            % Which is the first intersection point based off of vector direction
            % https://math.stackexchange.com/questions/311921/get-location-of-vector-circle-intersection
            a = (p1(1)-p0(1))^2 + (p1(2)-p0(2))^2;  % Quadratic "a"
            b = 2*(p1(1)-p0(1))*(p0(1)-center(1)) + ...
                2*(p1(2)-p0(2))*(p0(2)-center(2));                      % Quadratic "b"
            c = (p0(1)-center(1))^2 + (p0(2)-center(2))^2 - SOI_rad^2;  % Quadratic "c"
            t = (2*c)/(-b + sqrt(b^2 - 4*a*c));                         % Quadratic roots
    
            % Plug the root t back into the parametric form of the line
            inter_true(1) = (p1(1)-p0(1))*t + p0(1);     
            inter_true(2) = (p1(2)-p0(2))*t + p0(2);
            
            
            if any(inter_true == inter_1_pos)
                
                % Get the linear interpolated velocity values
                [~,~,ratio] = linearInterpolation(p0,p1,inter_1_pos(1));
                pv0 = [sc_vX(i2), sc_vY(i2)];      % Point 1 of the line
                pv1 = [sc_vX(i2+1), sc_vY(i2+1)];  % Point 2 of the line
                [vX_new,vY_new,~] = linearInterpolation(pv0,pv1,NaN,ratio);
                inter_1_vel = [vX_new,vY_new];  
                
                return_values = [inter_true;inter_1_vel];
                return
                
            elseif any(inter_true == inter_2_pos)
                
                % Get the linear interpolated velocity values
                [~,~,ratio] = linearInterpolation(p0,p1,inter_2_pos(1));
                pv0 = [sc_vX(i2), sc_vY(i2)];      % Point 1 of the line
                pv1 = [sc_vX(i2+1), sc_vY(i2+1)];  % Point 2 of the line
                [vX_new,vY_new,~] = linearInterpolation(pv0,pv1,NaN,ratio);
                inter_2_vel = [vX_new,vY_new];
                
                return_values = [inter_true;inter_2_vel];
                return
            
            else
                errorPathDisplay();
                fprintf(2,'Intersection points not matching up in SOIIntersectLocation.m\n')
                return
            end

            
        else
            continue

        end
        
    end
    
% end



end


