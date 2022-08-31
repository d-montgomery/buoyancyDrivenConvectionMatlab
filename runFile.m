%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Buoyancy Driven Convection in a Box:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code simulates buoyancy driven convection by placing a cold plate on
% top of a box filled with hot water. This is done by using Boussinesq's 
% approximation and solving the Navier-Stokes equations with
% the pressure from that approximation, and solving the unsteady
% advection-diffusion equation for the temperature.
%
% 2D Navier-Stokes
%       rho[ u_t + u u_x + v u_y] =  - p_x + mu ( u_xx + u_yy)
%       rho[ v_t + u v_x + v v_y] =  - p_y + mu ( v_xx + v_yy) - rho g
% for x in [0, Lxb], y in [0, Lyb] with Dirichlet BCs
%
% Boussinesq's approximation: 
%       rho g = rho_o [1 - B(T - Th)] g, where B = 207e-6 K^{-1}.
% Subbing into y-comp of Navier-Stokes (and doing some math) results in
% rho[ v_t + u v_x + v v_y] =  - (p + rho_o g y)_y + mu ( v_xx + v_yy) 
%                              + rho_o B (T - Th) g 
% 
% Unsteady Advection-Diffusion
%       T_t + u T_x + v T_y = alpha (T_xx + T_yy),
%   with Dirichlet conditions on top/bottom and zero flux on left/right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;

% Grid flag: 0 -> uniform grid
%            1 -> Chebyshev Points
% Plot flag: 1 -> Plot of the FV Grid (try Nx = Ny = 5 and Tf = dt), 
%            2 -> Surf during time stepping
% incr:      Determine how often plots are generated (if Plot flag = 2)
% Save flag: 0 -> don't save figures
%            1 -> save movie of streamlines and heatmap
    
gridFlag = 1;
plotFlag = 2; 
incr = 10;
saveFlag = 0;

% Domain = [0, Lx] x [0, Ly] ( meters )
Lx = 0.01; 
Ly = 0.01;
Nx = 200; % Number of Finite Volume cells in x-dir
Ny = Nx; % Number of Finite Volume cells in y-dir

% Temporal Parameters 
dt = 0.001; % Times step in seconds
Tf = 10; % Final time in seconds

% run main solver (! DON'T EDIT THESE 2 LINES !)
flags = [plotFlag, saveFlag, gridFlag, incr]; 
main_solver(Nx,Ny,Lx,Ly,Tf,dt,flags)