function [a,b,BC] = boundary_conditions(Ly)

% Get necessary parameters
[~,~,~,Tc,Th,~,~,~] = parameters(Ly);

% RHS BC's for x-equation (Navier-Stokes)
% u = g on each boundary (1 left, 2 right, 3 top, 4 bottom)
g1 = @(y,t) 0*y;
g2 = @(y,t) 0*y;
g3 = @(x,t) 0*x;
g4 = @(x,t) 0*x;
    
% RHS BC's for y-equation (Navier-Stokes)
% v = h on each boundary (1 left, 2 right, 3 top, 4 bottom)
h1 = @(y,t) 0*y;
h2 = @(y,t) 0*y;
h3 = @(x,t) 0*x;
h4 = @(x,t) 0*x;


% Set Boundary Conditions for Temp-equation ------------------------------
% General BC: a*T + b*grad(T) \cdot n = w 
% Note: a1 = 0, b1 = 1 provides T_x = w1 on left boundary
%       a2 = 0, b2 = 1 provides T_x = w2 on right boundary
%       a3 = 1, b3 = 0 provides T = w3 on top boundary
%       a4 = 1, b4 = 0 provides T = w4 on bottom boundary
a1 = 0; a2 = 0; a3 = 1; a4 = 1;
b1 = 1; b2 = 1; b3 = 0; b4 = 0;

% RHS BC's for T-equation
w1 = @(y,t) 0*y;
w2 = @(y,t) 0*y;
w3 = @(x,t) Th + 0*x; 
w4 = @(x,t) Tc + 0*x; 

% Create cell and vector for Boundary Conditions 
BC = {g1, g2, g3, g4;
      h1, h2, h3, h4;
      w1, w2, w3, w4};

a = [a1 a2 a3 a4];
b = [b1 b2 b3 b4];