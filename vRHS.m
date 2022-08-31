function[F] = vRHS(Nx,Ny,V0,NLV,NLV0,P,Qy,Qyo,dt,rho,mu,xu,yu,xv,yv,BC)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% V0 is velocity V from previous time step
% NLV and NLV0 are matrices of the Nonlinear terms for V
% P0 is pressure from the previous time step
% Qy and Qyo are the body forcing for the current and prev. step
% dt is time step
% rho is density
% mu is viscosity constant
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% BC is a cell arraw with RHS vectors for BC: v = h

h1 = BC{1};
h2 = BC{2};
h3 = BC{3};
h4 = BC{4};

F = zeros( (Ny+1)*(Nx+2),1 );

% lower/upper walls
for j = 2:Nx+1  
  F(j) = h3(j); % lower wall
  F((Nx+2)*Ny + j) = h4(j); % upper wall
end

% left/right boundaries
for i = 1:Ny+1
  F( (i-1)*(Nx+2)+1 ) = h1(i);   % left boundary
  F( i*(Nx+2) ) = h2(i);   %right boundary
end

% interior points
for i = 2:Ny
for j = 2:Nx+1
    
    % Cell height/width
    dx = xu(j) - xu(j-1);
    dy = yu(i+1) - yu(i);
    
    % Location of Nodes
    XE = xv(j+1);  XP = xv(j);  XW = xv(j-1);
    YN = yv(i+1);  YP = yv(i);  YS = yv(i-1);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -mu*( VP - VW )*dy/(XP-XW)/2;
    de = mu*( VE - VP )*dy/(XE-XP)/2;
    ds = -mu*( VP - VS )*dx/(YP-YS)/2;
    dn = mu*( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;
    
    % Pressure
    Pn = P(i,j-1); Ps = P(i-1,j-1);
    PR = dx * (Pn-Ps);
    
    F( (i-1)*(Nx+2) + j ) = rho*V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + (Qy(i,j) + Qyo(i,j))*dy*dx/2;
  
end
end






