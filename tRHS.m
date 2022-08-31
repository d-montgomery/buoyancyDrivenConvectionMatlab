function[F] = tRHS(Nx,Ny,T0,NL,NL0,Qs,Qso,xt,yt,xu,yv,dt,alpha,BC)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% T0 is temp from previous time step
% NL and NL0 are matrices of the Nonlinear terms for T
% Qs and Qso are the body forcing for the current and prev. step
% xt and yt is where T(x,y,t) is stored on grid
% xu is the x-comp of where u(x,y) is stored on grid
% yv is the y-comp of where v(x,y) is stored on grid
% dt is temporal step size
% alpha is the thermal diffusivity coefficient
% BC is a cell arraw with RHS vectors for BC: a T + b grad(T) \cdot n = w

% OUTPUT:
% F is a vector

w1 = BC{1};  w2 = BC{2};  w3 = BC{3};  w4 = BC{4};

% Initialize the RHS vector, F
F = zeros( (Ny+2)*(Nx+2),1 );

% left/right boundaries
for i = 1:Ny+2
  
  F( (i-1)*(Nx+2)+1 ) = w1(i);   % left boundary
  F( i*(Nx+2) ) = w2(i);   %right boundary

end

% lower/upper boundaries
for j = 2:Nx+1
  
  F(j) = w3(j); % lower boundary
  F((Nx+2)*(Ny+1) + j) = w4(j); % upper boundary
  
end

% interior points
for i = 2:Ny+1
for j = 2:Nx+1
  
    dx = xu(j) - xu(j-1);  
    dy = yv(i) - yv(i-1);
    
    % compute x and y locations of nodes in stencil
    xe = xt(j+1);
    xj = xt(j);
    xw = xt(j-1);
    yn = yt(i+1);
    yi = yt(i);
    ys = yt(i-1);
    
    % compute diffusion terms
    de = alpha*dy*( T0(i,j+1) - T0(i,j) )/(xe-xj);
    dn = alpha*dx*( T0(i+1,j) - T0(i,j) )/(yn-yi);
    dw = -alpha*dy*( T0(i,j) - T0(i,j-1) )/(xj-xw);
    ds = -alpha*dx*( T0(i,j) - T0(i-1,j) )/(yi-ys);
    CN = 0.5*(de+dn+dw+ds);
    
    F( (i-1)*(Nx+2) + j ) = T0(i,j)*dx*dy/dt -1.5*NL(i,j) +0.5*NL0(i,j)...
                            + CN + 0.5*(Qs(i,j) + Qso(i,j))*dy*dx;
    
end
end