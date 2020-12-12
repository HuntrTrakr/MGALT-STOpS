%House Cleaning
close all; clear; clc

%Add other folders
addpath Error_Messages\                 ...
        Optimization_Results\           ...
        Other_Functions\                ...
        Solvers\                     	...
        Solvers\Cost_Functions\         ...
        Solvers\Direct\               	...
        Solvers\Indirect\           	...
        Solvers\Transfer_Conditions\	...
        Thrust_Profiles
    

    
%% Current Verification Results

load DEMO_results.mat;

% load your_filename_here.mat;



%% Display Results

displayResults(BOD,CONST,OPT,VAR,run_time,eval_info)


