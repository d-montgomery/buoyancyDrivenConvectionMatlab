# Buoyancy Driven Convection 
This code simulates buoyancy driven convection using Boussinesq's approximation.  The simulation is setup so that a box of water at rest has an initial temp of Th, where the bottom of the box is fixed at the temp Th. Gravity is the only force acting on the liquid and g is defined to act only in the negative y direction.  At time t = 0, a cold plate is placed on top of the box with temp Tc. If the temperature between the two plates (Th and Tc) is large enough, the cold will be advected and mixed.  If the temperature between the two plates is small, the cold will be slowly diffuse down. To that end, this code has been setup to simulate the case where the temperature difference between Tc and Th is large enough to create convective motion.

Running this code: 
To run this code, simply open runFile.m in matlab and run. 

Editing code:
1) runFile.m allows you to edit the spatial and temporal paramters.  Additionally, there are different "flags" (lines 36-39) that allow you to determine what type of grid you'd like to use, what type of plot you'd like to create, how often the plots should be made, and whether or not you'd like to create a movie of the simulation.
2) paramters.m allows you to edit all of the physical parameters.  If you'd like to change the type of fluid, simply edit lines 5-8.  You can also change the Rayleigh number and temperatures of the bottom and top plates in this file.

If you have any questions feel free to email me at d.montgomery016@gmail.com

 
