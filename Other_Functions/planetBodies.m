function [bodies_R,bodies_V,bodies_JD] = planetBodies(BOD)
% FORM: [bodies_R,bodies_V,bodies_JD] = planetBodies(BOD)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Return R, V, and JD arrays for bodies included in the problem for 
% |     the valid time span
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -bodies_R     	(3*Nbodies,n)	[float]         [AU]
% |         Heliocentric Radius components for all planets
% |     -bodies_V       (3*Nbodies,n)	[float]         [AU/TU]
% |         Heliocentric Velocity components for all planets
% |     -bodies_JD   	(Nbodies,n)     [float]         [JD]
% |         Julian day for the departure body and all transfer bodies
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Initial

bodies = BOD.bodies;
tspan = BOD.tspan;


% Calculate Start/End year based off tspan
[year_span(1),~,~,~,~,~,~,~] = julian2greg(tspan(1));   % Bug fix...kept making a JD for May 6 2021 as 2022
[year_span(2),~,~,~,~,~,~,~] = julian2greg(tspan(2));

% Preallocate
Data = cell(length(bodies),(year_span(2)-year_span(1)+1));
total_indicies = 0;
array_indicies = [0,0];

% If within the same year
if year_span(1) == year_span(2)
    year_span(2) = [];
end

% Add the years between the start and ending year
if length(year_span) > 1    % If larger than 1 year
    
    year = zeros(1,(year_span(2)-year_span(1)+1));
    
    for i1 = 1:length(year)
        year(i1) = year_span(1)+i1-1;
    end
    
else
    year = year_span;
end

% Load Appropriate Data
for i2 = 1:length(year)
    
    for i3 = 1:length(bodies)
        
        string = [bodies{i3},'_',num2str(year(i2))];
        Data{i3,i2} = load(string);
        % Add something for when data is not available (ode45)
    end
    
end

% Make a array of the matrix indicies for output data
for i4 = 1:length(year)
    
    [~,temp_array] = size(Data{1,i4}.R);
    array_indicies(i4,:) = [array_indicies(end,end)+1,...
        array_indicies(end,end)+temp_array];
    total_indicies = total_indicies + temp_array;
    
end

% Pre allocate
bodies_R = zeros(length(bodies)*3,total_indicies);
bodies_V = zeros(length(bodies)*3,total_indicies);
bodies_JD = zeros(1,total_indicies);

% Put the data from the extracted files into the return values
for i5 = 1:length(bodies)
    
    for i6 = 1:length(year)
        
        rows = i5*3-2:i5*3;
        cols = array_indicies(i6,1):1:array_indicies(i6,end);
        bodies_R(rows,cols) = Data{i5,i6}.R;
        bodies_V(rows,cols) = Data{i5,i6}.V;
        
        if i5 == 1
            bodies_JD(i5,cols) = Data{i5,i6}.JD;
        end
        
        clear cols
    end
    
end



end


