function [num] = randomNum(a,b,dec_int)
% FORM: [num] = randomNum(a,b,dec_int)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function generates a random number between the given 
% |     boundaries. The output can be either an integer value, or a 
% |     floating point decimal value. If the choice for the output is 
% |     'int', the boundaries "a" and "b" must be integers; otherwise, 
% |     they will be rounded to the nearest integer.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -a                  (1,1)     	[float]         [unitless]
% |         Lower boundary
% |     -b                  (1,1)     	[float]         [unitless]
% |         Upper boundary
% |     -dec_int          	(1,n)      	[string]    	[unitless]
% |         Choice between integer or decimal value for chosen random 
% |         number. Choices: 'int' or 'dec'
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -num                (1,1)       [float][int]    [unitless]
% |         Randomly chosen number between (a,b)
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Make random number

switch dec_int
    
    case {'int'}
        
        a = round(a) - 0.5;
        b = round(b) + 0.5;
        num = round(a+(b-a)*rand); 
        
    case {'dec'}
        
        num = a+(b-a)*rand;
        
    otherwise
        
        errorPathDisplay();
        fprintf(2,'Incorrect random choice selected.\n\n')
        fprintf(2,'Valid options are:\n')
        fprintf(2,'"int"\n')
        fprintf(2,'"dec"\n')
        return
        
end



end


