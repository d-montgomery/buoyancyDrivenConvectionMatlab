function [Qx,Qy,Qt] = body_forcing(XU,XT,T0v,Ly)

% Get Necessary Parameters
[~,rho,~,~,Th,~,B,g] = parameters(Ly);

% Set Body Force Terms Qx, Qy, Qt
Qx = 0 * XU;
Qy = rho * B * g * (T0v - Th);
Qt = 0 * XT;
          