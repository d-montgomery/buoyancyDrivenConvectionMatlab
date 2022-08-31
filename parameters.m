function [mu,rho,alpha,Tc,Th,STT,B,g] = parameters(Ly)
% Input: Ly is the height of the box (needed for Th)

% Physical parameters for the fluid (water)
mu = 10^(-3); % kg/(ms) (diffusion coefficient)
rho = 998.2; % kg/(m^3) (density)
alpha = 0.143e-6; % m^2/s (thermal diffusivity)
B = 207e-6; % K^-1 (thermal expansion)
g = 9.8; %ms^-2 (Gravity acting on fluid)

% Parameters for temp
Tc = 293.15; % K (temp of cold plate on top)
Ra = 4e6; % Rayleigh number (fixed so Delta T can be solved for)

DT = Ra*alpha*mu/(B*g*rho*Ly^3); % proper Delta T for buoyancy driven conv.
Th = Tc + DT; % K (temp of bottom plate) 
STT = 273.15; % Standard temp in K

