function [JD] = getJulianDate(date)
% FORM: [JD] = getJulianDate(date)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -Get the Julian Date without a toolbox
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -date               (1,1)       [string]        [date]
% |         String of date to convert in YYYY/MM/DD
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -JD                 (1,1)       [float]         [JD]
% |         The Julian Date
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     https://www.mathworks.com/matlabcentral/answers/1303-convert-julian-date-to-calendar-days
% |
% |-----------------------------------------------------------------------



%% Get JD

JD = datenum(date) + 1721058.5;



end


