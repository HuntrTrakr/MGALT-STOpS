MGALT-STOpS 1.2
===============

MGALT-STOpS is a software suite which allows for optimization of low-thrust spacecraft trajectories between multiple planets.
The Spacecraft Trajectory Optimization Suite (STOpS) uses Multiple Gravity-Assist Low-Thrust (MGALT) trajectories paired with the island model paradigm to accomplish this goal. 
The island model utilizes four different global search algorithms: a Genetic Algorithm, Differential Evolution, Particle Swarm Optimization, and Monotonic Basin Hopping.
Each island runs either a direct or indirect optimization method to solve for the most optimal trajectory between multiple planets.


## Examples

Included in this download is a demo program, titled "MGALT_STOpS_DEMO", which provides a step by step walkthrough of MGALT-STOpS and explains everything along the way.
It is recommended to become familiar with the demo script before venturing into customization with the main script.


## Framework

All of the files were built and tested with MATLAB 2019b [MATLAB](https://www.mathworks.com/products/matlab.html).
Future versions of MATLAB may not be supported, and there are no plans for continued development of MGALT-STOpS on future versions of MATLAB.

No toolboxes are necessary to run MGALT-STOpS; however, the "Parallel Computing in Optimization Toolbox" can be used to greatly speed up execution time of MGALT-STOpS.
There is an option to enable/disable parallel computing at the beginning of the scripts.


## This Version

Version 1.2 fixes a bug with the implementation of parallel computing in the MBH secondary search, updates default weighting parameters, and makes the arrows on the direct method plot more consistent.


## More Information

This code was inherited/developed/updated by Malloy to support his thesis, "Spacecraft Trajectory Optimization Suite (STOpS): Design and Optimization of Multiple Gravity-Assist Low-Thrust (MGALT) Trajectories Using Modern Optimization Techniques"; a requirement for a MS in Aerospace Engineering at California Polytechnic State University, San Luis Obispo.

The development of MGALT-STOpS was a continuation of prior work on STOpS by Shane P. Sheehan and Timothy J. Fitzgerald.
The supporting thesis for MGALT-STOpS, Low-Thrust STOpS, and STOpS can be found below:
* ["Spacecraft Trajectory Optimization Suite (STOpS): Design and Optimization of Multiple Gravity-Assist Low-Thrust (MGALT) Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/2247/) by MALLOY (2020)
* ["Spacecraft Trajectory Optimization Suite (STOPS): Optimization of Low-Thrust Interplanetary Spacecraft Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/1901/) by SHEEHAN (2017)
* ["Spacecraft Trajectory Optimization Suite (STOpS): Optimization of Multiple Gravity Assist Spacecraft Trajectories Using Modern Optimization Techniques"](https://digitalcommons.calpoly.edu/theses/1503/) by FITZGERALD (2015)

The author can be contacted on Discord at "HuntrTrakr#8835" if there are any questions.
