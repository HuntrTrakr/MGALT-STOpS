function [] = errorPathDisplay()
% FORM: [] = errorPathDisplay()
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Displays a simple message telling the user that an error was
% |     detected. The function then displays the function stack in the 
% |     command window for easy debugging
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
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



%% Info

% Get the location and line numbers of the function stack
error_msg = dbstack();

% Print the warning message!
fprintf(2,'\n\n\n---------------------------------------\n')
fprintf(2,'WARNING! Error detected.\n')
fprintf(2,'Error stack:\n\n')
fprintf(2,'Line:      File:\n')

% Print off the function call stack
for i1 = 2:size(error_msg,1)
    
    fprintf(2,'%3.0f   in   %s\n',error_msg(i1,1).line,error_msg(i1,1).name)
    
end

fprintf(2,'---------------------------------------\n\n\n\n')



end


