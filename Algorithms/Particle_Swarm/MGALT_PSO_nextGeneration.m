function [bee] = MGALT_PSO_nextGeneration(BOD,OPT,OPT_algo,VAR,bee)
% FORM: [bee] = MGALT_PSO_nextGeneration(BOD,OPT,OPT_algo,VAR,bee)
%
% |-----------------------------------------------------------------------
% |
% | NOTES:
% |     -This function creates the next generation of PSO popn from the 
% |     old member population and slight perturbations to the old members.
% |
% |-----------------------------------------------------------------------
% |
% | INPUTS:
% |     -BOD                (1,1)       [struct]        [unitless]
% |         A struct containing information pertaining to the planetary
% |         bodies. Contains list of bodies, launch windows and ToF, and 
% |         planetary R/V/JD vectors. This struct has dynamic fields and 
% |         will adapt to contain only the necesary information
% |     -OPT                (1,1)       [struct]        [unitless]
% |         A struct containing constants user options. Contains the save 
% |         folder, ToF values, and more structs containing informaiton 
% |         for the island model, cost parameters, weighting parameters, 
% |         and all of the islands used in the optimization process
% |     -OPT_algo           (1,1)       [struct]        [unitless]
% |         PSO option parameters. For a full explination of these 
% |         parameters, see 
% |         "Algorithms/Algorithm_Parameters/parametersPSO.m"
% |     -VAR                (1,1)       [struct]        [unitless]
% |         A struct containing the variable limits
% |     -bee                (Npop,1)	[struct]       	[unitless]
% |         Contains cost, speed position, ect of every bee in the 
% |         population for the old generation
% |
% |-----------------------------------------------------------------------
% |
% | OUTPUTS:
% |     -bee                (Npop,1)	[struct]       	[unitless]
% |         Contains cost, speed position, ect of every bee in the 
% |         population for the new generation
% |
% |-----------------------------------------------------------------------
% |
% | MISC:
% |
% |-----------------------------------------------------------------------



%% Setup

% Assign Sensible Names to Variables
Vmax = OPT_algo.vmax*(VAR.high-VAR.low);
K = OPT_algo.K;
c1 = OPT_algo.cl;
cmax = OPT_algo.cmax;
Npop = length(bee);
Nvar = length(bee(1).pos);
f_inform = zeros(K,1);



%% Propegate Bees

for b = 1:length(bee)
   informant = zeros(K,1);
   
   % Determine Informants
   for I = 1:K
      while informant(I) == 0 || informant(I) == b || informant(I) == any(informant(1:(I-1)))
         informant(I) = randomNum(1,Npop,'int');
      end
      f_inform(I) = bee(informant(I)).f;
   end
   
   % Determine Best Informant
   [ ~, ind ] = min(f_inform);
   if f_inform(ind) < bee(b).f_g
      bee(b).g   = bee(informant(ind)).pos;
      bee(b).f_g = f_inform(ind);
   end
   
   % Move Member & Re-Calculate Velocity
   bee(b).vel = c1*bee(b).vel + rand*cmax*(bee(b).p-bee(b).pos) + rand*cmax*(bee(b).g-bee(b).pos);
   bee(b).pos = bee(b).pos + bee(b).vel;

   % Ensure Bounds Not Exceeded
   for i = 1:Nvar
       
      if bee(b).vel(i) >  Vmax(i)
          bee(b).vel(i) =  Vmax(i);
      end
      
      if bee(b).vel(i) < -Vmax(i)
          bee(b).vel(i) = -Vmax(i);
      end
      
      if bee(b).pos(i) >  VAR.high(i)
          bee(b).pos(i) = VAR.high(i); 
          bee(b).vel(i) = -bee(b).vel(i);
      end
      
      if bee(b).pos(i) <  VAR.low(i)
          bee(b).pos(i) = VAR.low(i); 
          bee(b).vel(i) = -bee(b).vel(i); 
      end
      
      if VAR.bin(i)
          bee(b).pos(i) = round(bee(b).pos(i));
      end

   end
   
   % Ensure ToF is not broken
   bee(b).pos = MGALT_fixToF(BOD,OPT,VAR,bee(b).pos);
   
end
   


end


