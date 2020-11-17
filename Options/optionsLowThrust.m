function [opt_thrust] = optionsLowThrust(OPT,selection)
% FORM: [opt_thrust] = optionsLowThrust(OPT,selection)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function call outputs the necessary struct object to fully
% |     define a Differential Evolution Island. A switch/case allows the 
% |     user to add in commonly used cases to have it as an easy return, 
% |     and more cases can easily be added. The function also allows for 
% |     custom inputs, which is then exported to be in the correct form 
% |     for the Island struct object.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -selection          (1,n)       [string]        [unitless]
% |         A used defined string which is used in a switch/case format
% |         to return predefined options
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -opt_thrust         (1,n)       [assorted]      [unitless]
% |         The return object with all the associated parameters
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |     The parameters are defined as follows:
% |
% |
% |     For ALL Methods:
% |
% |     -tt_end:            (1,1)       [float]         [days]
% |         Target end time. This is approximately the desired end time 
% |         for the trajectory in days. it is important for the INDIRECT 
% |         method this is short enough or a good solutions will not be 
% |         found.
% |     -time:              (1,2)       [int]           [unitless]
% |         Margin on end time. Basically the upper and lower bounds for 
% |         the transfer end time with the first value being the lower 
% |         bound and the second value being the upper bound. The end 
% |         time of the trajectory is allowed to vary from tt_end-time(1) 
% |         to tt_end+time(2).
% |     -thrust_method:     (1,1)       [string]        [unitless]
% |         Thrust method. Options are: 'constant', 'variable', and 
% |         'equation'. The INDIRECT method can only handle 'constant' 
% |         and 'equation'. All entries are spelling and case sensitive.
% |     -thrust:            (1,1)(1,2) [float][string]  [unitless]
% |         Thrust for the trajectory. If thrust method is 'constant' 
% |         this is a (1,1) [#] for the thrust value in Newtons. If the 
% |         thrust method is 'variable_thrust' then this is a (1,2) [#] 
% |         for upper and lower bounds on the thrust in Newtons. If the 
% |         thrust method is 'equation_thrust' this is a function handle 
% |         to calculate the thrust; inputs to a thrust function must be 
% |         radius (AU) and Isp (s) in that order. Function handle is 
% |         spelling and case sensitive. 
% |     -m0:                (1,1)       [float]         [kg]
% |         Inital wet mass of spacecraft in kilograms.
% |     -mdot_method:       (1,1)       [string]        [unitless]
% |         Method for calculating mass flow rate. Options are: 
% |         'constant' and 'equation'. All options are spelling and case 
% |         sensitive.
% |     -mdot:              (1,1)       [float][string]	[kg/s]
% |         Mass flow rate for the trajectory. If mdot method is 
% |         'constant' this is the mass flow rate in kg/s. If the mdot 
% |         method is 'equation_mdot' this is the function handle to 
% |         calculate the mdot; inputs to a mdot equation must be 
% |         thrust (N) and Isp (s) in that order. Function handle is 
% |         spelling and case sensitive.
% |     -Isp:               (1,1)       [float]         [seconds]
% |         The Isp of the propulsion method in seconds.
% |
% |
% |
% |     For the DIRECT Method
% |
% |     -Nseg:              (1,1)       [int]           [unitless]
% |         Number of segments that the trajectory is divided into. The 
% |         recommended starting values is ten. This is only valid for 
% |         the DIRECT Segmented method.
% |     -orbit_check:       (1,1)       [string]        [unitless]
% |         With MALLOY'S version of MGALT STOpS, having the sting on
% |         anything except for 'off' will result in an error!
% |         Whether or not the number of heliocentric revolutions before 
% |         converging on the target is controlled. Entering 'on' allows 
% |         the user to specify a number of heliocentric orbits for the 
% |         spacecraft to complete before attempting to converge on the 
% |         destination. Entering 'off' means the suite will attempt to 
% |         find the best number of orbits on its own.
% |     -orbits:            (1,1)       [int]           [unitless]
% |         The number of heliocentric orbits required before 
% |         attempting to converge on the target planet/orbit if 
% |         orbit_check = 'on'. Entering 1 means the spacecraft will 
% |         complete 0-360 degrees of heliocentric orbit before 
% |         converging. Entering 2 is 360-720 degrees of heliocentric 
% |         orbit, ect.
% |
% |-----------------------------------------------------------------------



%% Selection

opt_thrust = struct;

switch selection
    
    case {'Custom-1'}
        
        opt_thrust.Nseg = 10;
        opt_thrust.orbit_check = 'off';
        opt_thrust.orbits = 0;
        opt_thrust.tt_end = OPT.tof_total;
        opt_thrust.time = OPT.tof_margin;
        opt_thrust.m0 = 15000;
        opt_thrust.thrust_method = 'constant';
        opt_thrust.thrust = 1.5;
        opt_thrust.mdot_method = 'constant';
        opt_thrust.mdot = 20e-6;
        opt_thrust.Isp = 5000;
        opt_thrust.n_available = 2;
        opt_thrust.duty_cycle_type = 'calculated';
        opt_thrust.duty_cycle = 1;
        opt_thrust.launch_dV_rad = 0;
        opt_thrust.launch_dV_tan = 0;
    
    case {'Rauwolf'}
        
        opt_thrust.Nseg = 10;
        opt_thrust.orbit_check = 'off';
        opt_thrust.orbits = 0;
        opt_thrust.tt_end = OPT.tof_total;
        opt_thrust.time = OPT.tof_margin;
        opt_thrust.m0 = 4545.5;
        opt_thrust.thrust_method = 'constant';
        opt_thrust.thrust = 3.787;
        opt_thrust.mdot_method = 'equation';
        opt_thrust.mdot = 'SolarSail_mdot';
        opt_thrust.Isp = 3000;
        opt_thrust.n_available = 1;
        opt_thrust.duty_cycle_type = 'calculated';
        opt_thrust.duty_cycle = 1;
        
    case {'Yam-STOUR'}
        
        opt_thrust.Nseg = 30;
        opt_thrust.orbit_check = 'off';
        opt_thrust.orbits = 0;
        opt_thrust.tt_end = OPT.tof_total;
        opt_thrust.time = OPT.tof_margin;
        opt_thrust.m0 = 20000;
        opt_thrust.thrust_method = 'constant';
        opt_thrust.thrust = 2.26;
        opt_thrust.mdot_method = 'constant';
        opt_thrust.mdot = 38.4e-6;
        opt_thrust.Isp = 6000;
        opt_thrust.n_available = 5;
        opt_thrust.duty_cycle_type = 'calculated';
        opt_thrust.duty_cycle = 1;
        opt_thrust.launch_dV_rad = 0;
        opt_thrust.launch_dV_tan = 0.2;
        
    case {'Yam-GALLOP'}
        
        opt_thrust.Nseg = 30;
        opt_thrust.orbit_check = 'off';
        opt_thrust.orbits = 0;
        opt_thrust.tt_end = OPT.tof_total;
        opt_thrust.time = OPT.tof_margin;
        opt_thrust.m0 = 20000;
        opt_thrust.thrust_method = 'constant';
        opt_thrust.thrust = 2.26;
        opt_thrust.mdot_method = 'constant';
        opt_thrust.mdot = 38.4e-6;
        opt_thrust.Isp = 6000;
        opt_thrust.n_available = 5;
        opt_thrust.duty_cycle_type = 'calculated';
        opt_thrust.duty_cycle = 1;
        opt_thrust.launch_dV_rad = 0;
        opt_thrust.launch_dV_tan = 0;
        
    otherwise                   	% Error
        
        errorPathDisplay();
        return
        
end



end


