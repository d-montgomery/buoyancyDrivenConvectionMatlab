function[F] = uRHS(Nx,Ny,U0,NLU,NLU0,P0,Qx,Qxo,dt,rho,mu,xu,yu,xv,yv,BC)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% U0 is velocity U from previous time step
% NLU and NLU0 are matrices of the Nonlinear terms for U
% P0 is pressure from the previous time step
% Qx and Qxo are the body forcing for the current and prev. step
% dt is time step
% rho is density
% mu is viscosity constant
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% BC is a cell arraw with RHS vectors for BC: u = g

% OUTPUT:
% F is a vector

g1 = BC{1}; 
g2 = BC{2};
g3 = BC{3};
g4 = BC{4};

F = zeros( (Ny+2)*(Nx+1),1 );

% lower/upper walls
for j = 2:Nx
  F(j) = g3(j); % Lower wall
  F((Nx+1)*(Ny+1)+j) = g4(j); % upper wall
end

% left/right boundaries
for i = 1:Ny+2
  F( (i-1)*(Nx+1)+1 ) = g1(i);   % left boundary
  F( i*(Nx+1) ) = g2(i);        %right boundary
end

% interior points
for i = 2:Ny+1
for j = 2:Nx
    
    % Cell height/width
    dx = xv(j+1) - xv(j);
    dy = yv(i) - yv(i-1);
    
    % Location of Nodes
    XE = xu(j+1);  XP = xu(j);  XW = xu(j-1);
    YN = yu(i+1);  YP = yu(i);  YS = yu(i-1);
    
    % Velocities           
    UP = U0(i,j);
    UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


    % Diffusion Terms
    dw = -mu*(UP - UW)*dy/(XP-XW)/2;
    de = mu*(UE - UP)*dy/(XE-XP)/2;
    ds = -mu*(UP - US)*dx/(YP-YS)/2;
    dn = mu*(UN - UP)*dx/(YN-YP)/2;
    D = de + dn + dw + ds;

    % Pressure Terms
    Pe = P0(i-1,j); Pw = P0(i-1,j-1);
    PR = dy * (Pe-Pw);

    F( (i-1)*(Nx+1) + j ) = rho*U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                              0.5*NLU0(i,j) + D - PR + ...
                             + (Qx(i,j) + Qxo(i,j))*dx*dy/2;
  
end  
end






