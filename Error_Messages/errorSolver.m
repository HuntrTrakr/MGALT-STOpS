function [] = errorSolver()
% FORM: [] = errorSolver()
%
% |-----------------------------------------------------------------------
% | NOTES:
% |     -Displays a simple message telling the user that an invalid
% |     solver was selected and displays the valid choices in the 
% |     command window
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

fprintf(2,'Incorrect solver selected.\n\n')
fprintf(2,'Valid options are:\n')
fprintf(2,'LT_DIR_FSM_2D\n')
fprintf(2,'LT_IN_FSM_2D\n')
fprintf(2,'MGALT_DIR_FBSM_2D\n')
fprintf(2,'MGALT_IN_FBSM_2D\n\n')
% fprintf(2,'"NEW SOLVER"\n')



end


