# Four-state Monte Carlo model of the microtubule dynamics

Matlab2017b scripts for the four-state Monte Carlo model of microtubule dynamics [ref: Alexandrova et al., 2022. A theory of tip structure-dependent microtubule catastrophes and damage-induced microtubule rescues]: 
-	_MT_4state_simulation.m_ is the main script used to run the simulations;
-	_Visualization_MTline.m_ is the script used to visualize microtubules in the simulations.

## Configuring and running a simulation
Open _MT_4state_simulation.m_ file and provide input parameters by changing the lines 25-37 in MT_4state_simulation.m. The default parameters are set according to the Table S1 (high activation barrier parameter set) of the paper (Alexandrova et al., 2022).  Run the script in Matlab to launch the simulation (1 minute of the simulated time is typically computed within 0.5-1 minutes on Intel Core i7 2.8 GHz processor). 

## Outputs
The program output includes plots of microtubule length, curved protofilament length, tip raggedness, GTP content and individual PFs length changes in the simulated trajectory. Also, in the command window one can find the average microtubule velocity (nm/s), average length of curved protofilament parts (nm), average raggedness of the microtubule (tip) over simulation time and the sum of GTP-subunits at the last iteration.

## Vizualization of one simulation frame
To visualize a microtubule at a given time point, run _Visualization_MTline.m_ file in Matlab at any moment during the simulation or in the end of the simulation. You can change the coordinates (in tubulin dimers) of the microtubule segment to be visualized, in the lines from 4 to 6. By default the program plots 30 terminal microtubule dimer layers. In line 7, one can change the projection angle of the microtubule (angle of rotation around the microtubule axis). The GTP-state is shown in red and GDP-state is green.
