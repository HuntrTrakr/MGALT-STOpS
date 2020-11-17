function [year,month,day,hour,min,sec,dayweek,dategreg] = julian2greg(JD)
% FORM: [year,month,day,hour,min,sec,dayweek,dategreg] = julian2greg(JD)
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -This function converts the Julian dates to Gregorian dates
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -JD                 (1,1)       [float]         [JD]
% |         The Julian Date to transform
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -year               (1,1)       [int]           [year]
% |         The year of the JD
% |     -month          	(1,1)       [int]           [month]
% |         The month of the JD
% |     -day                (1,1)       [int]           [day]
% |         The day of the JD
% |     -hour               (1,1)       [int]           [h]
% |         The hour of the JD in UT
% |     -min                (1,1)       [int]           [min]
% |         The min of the JD in UT
% |     -sec                (1,1)       [int]           [sec]
% |         The sec of the JD in UT
% |     -dayweek            (1,1)       [string]        [day]
% |         The day of the week of the JD
% |     -dategreg         	(1,6)       [int]           [date]
% |         The gergorian date of the JD
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     -Example
% |         [a,b,c,d,e,f,g,h] = julian2greg(2453887.60481)
% |         a = 2006      
% |         b = 6
% |         c = 1
% |         d = 2
% |         e = 30
% |         f = 56
% |         g = Thursday
% |         h = 1     6     2006     2     30     56
% |
% |     -References
% |         Astronomical Applications Department. "Julian Date Converter". 
% |             From U.S. Naval Observatory.
% |             http://aa.usno.navy.mil/data/docs/JulianDate.html
% |         Duffett-Smith, P. (1992).  
% |             Practical Astronomy with Your Calculator.
% |             Cambridge University Press, England:  pp. 8,9.
% |
% |     Gabriel Ruiz Mtz.
% |     Jun-2006
% |
% |-----------------------------------------------------------------------



%% Solve

narginchk(1,1)

I = floor( JD + 0.5);
Fr = abs( I - ( JD + 0.5) );	 

if I >= 2299160 
     A = floor( ( I- 1867216.25 ) / 36524.25 );
     a4 = floor( A / 4 );
     B = I + 1 + A - a4;
else
     B = I;
end 

C = B + 1524;
D = floor( ( C - 122.1 ) / 365.25 );
E = floor( 365.25 * D );
G = floor( ( C - E ) / 30.6001 );
day = floor( C - E + Fr - floor( 30.6001 * G ) );

if G <= 13.5 
    month = G - 1;
else
    month = G - 13;
end

if month > 2.5
    year = D - 4716;
else
    year = D - 4715;
end

hour = floor( Fr * 24 );
min = floor( abs( hour -( Fr * 24 ) ) * 60 );
minufrac = ( abs( hour - ( Fr * 24 ) ) * 60 ); 
sec = ceil( abs( min - minufrac ) * 60);
AA = ( JD + 1.5 ) / 7;
nd = floor( (abs( floor(AA) - AA ) ) * 7 );
dayweek = {'Sunday' 'Monday' 'Tuesday' 'Wednesday' 'Thursday' 'Friday' 'Saturday'};
dayweek = dayweek{ nd+1};
format('long', 'g');
dategreg = [ day month year hour min sec ];
	   

