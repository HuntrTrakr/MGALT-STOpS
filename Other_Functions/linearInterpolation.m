function [x,y,ratio_out] = linearInterpolation(p0,p1,x,ratio_in)
% FORM: [x,y,ratio_out] = linearInterpolation(p0,p1,x,ratio_in)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Function designed to find y given two line segments and a 
% |     desired x point somwhere between the lines. The funciton also 
% |     finds the ration of the distances between the initial point and 
% |     the desired point and will apply that to two other points to get 
% |     the linear interpolation from an unknown x value
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -p0                 (1,2)       [float]         [unitless]
% |         The first point with x and y locations
% |     -p1                 (1,2)       [float]         [unitless]
% |         The second point with x and y locations
% |     -x                  (1,1)       [float]         [unitless]
% |         The desired x value of a location somewhere between p0 and p1
% |     -ratio_in        	(1,1)       [float]       	[unitless]
% |         The ratio input as a previous calculation of the norm of
% |         magnitude of p0->x / p0->p1
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -x                  (1,1)       [float]         [unitless]
% |         The desired x value of a location somewhere between p0 and p1.
% |         This could be the input, or the unknown x value calculated by 
% |         the ratio_in
% |     -y                  (1,1)       [float]         [unitless]
% |         The desired y value of a point x
% |     -ratio_out        	(1,1)       [float]       	[unitless]
% |         The ratio as a  calculation of the norm of magnitude of 
% |         p0->x / p0->p1
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -Linear interpolation formula
% |         y = y0+(x-x0)*((y1-y0)/(x1-x0));
% |
% |-----------------------------------------------------------------------



%% Interpolate

% Formula for linear interpolation
y = p0(2) + (x-p0(1))*((p1(2)-p0(2))/(p1(1)-p0(1)));


% Get the magnitude of the start->end
mag_p0_p1 = sqrt( (p1(1)-p0(1))^2 + (p1(2)-p0(2))^2 );

% Get the magnitude from start->x
mag_p0_xy = sqrt( (p0(1)-x)^2 + (p0(2)-y)^2 );
mag_p0_xy_norm = mag_p0_xy/mag_p0_p1;

% Ratio of the magnitudes normalized to the start->end
ratio_out = mag_p0_xy_norm;


% If the ratio is an input
if exist('ratio_in','var')

    % Get the vector between start->end
    vec = [p1(1)-p0(1),p1(2)-p0(2)];
    
    % Applt the ratio to it
    vec = vec*ratio_in;

    % Get the new x from the ratio
    x = p0(1) + vec(1);
    
    % Linear interpolation for the new y
    y = p0(2) + (x-p0(1))*((p1(2)-p0(2))/(p1(1)-p0(1)));
    
end

end


