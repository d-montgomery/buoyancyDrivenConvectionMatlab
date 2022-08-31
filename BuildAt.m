function[A] = BuildAt(Nx,Ny,xt,yt,xu,yv,alpha,dt,a,b)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% xt and yt is where T(x,y,t) is stored on grid
% xu is the x-comp of where u(x,y) is stored on grid
% yv is the y-comp of where v(x,y) is stored on grid
% alpha is the thermal diffusivity coefficient
% dt is temporal step size
% a and b are vectors for the general BC: a T + b grad(T) \cdot n = w

% OUTPUT:
% A is a coefficient matrix

% Get the BC coefficients
a1 = a(1); a2 = a(2); a3 = a(3); a4 = a(4);
b1 = b(1); b2 = b(2); b3 = b(3); b4 = b(4);

% initialize the counter used in building the sparse matrix
ctr=0;

% initialize row, col and B. S is number of entries
S = 4*(Ny+2 - 1) + 4*(Nx+1-2) + 5*(Nx+1-2)*(Ny+1-2);
row = zeros(1,S);
col = row;
A = row;

% left boundary condition
dx = xt(2)-xt(1);
for i = 1:Ny+2
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2) + 1;
  col(ctr) = row(ctr);
  A(ctr)   = a1 - b1/dx;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2) + 1;
  col(ctr) = row(ctr)+1;
  A(ctr)   = b1/dx;
  
end

% right boundary condition
dx  = xt(Nx+2) - xt(Nx+1);
for i = 1:Ny+2
     
    ctr = ctr+1;
    row(ctr) = i*(Nx+2);
    col(ctr) = row(ctr);
    A(ctr)   = a2 + b2/dx;
    
    ctr = ctr+1;
    row(ctr) = i*(Nx+2);
    col(ctr) = row(ctr) -1;
    A(ctr)   = -b2/dx;
end

%lower wall boundary conditions
dy = yt(2)-yt(1);
for j=2:Nx+1
  
  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = j;
  A(ctr)   = a3-b3/dy;

  ctr = ctr+1;
  row(ctr) = j;
  col(ctr) = j + (Nx+2);
  A(ctr)   = b3/dy;
  
end

%upper wall boundary conditions
dy = yt(Ny+2)-yt(Ny+1);
for j=2:Nx+1
  
  ctr = ctr+1;
  row(ctr) = (Ny+1)*(Nx+2) + j;
  col(ctr) = row(ctr);
  A(ctr)   = a4 + b4/dy;
  
  ctr = ctr+1;
  row(ctr) = (Ny+1)*(Nx+2) + j;
  col(ctr) = row(ctr) - (Nx+2);
  A(ctr)   = -b4/dy;
    
end


% Internal Nodes
for j = 2:Nx+1
for i = 2:Ny+1   

  % compute size of cell faces
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
  dw = -alpha*dy/(xj-xw);
  de = -alpha*dy/(xe-xj);
  ds = -alpha*dx/(yi-ys);
  dn = -alpha*dx/(yn-yi);
  dp = -(dw+de+dn+ds);
  
  % apply coefficients
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j;
  A(ctr)   = dx*dy/dt + 0.5*dp;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j-1;
  A(ctr)   = 0.5*dw;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-1)*(Nx+2)+j+1;
  A(ctr)   = 0.5*de;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i)*(Nx+2)+j;
  A(ctr)   = 0.5*dn;
  
  ctr = ctr+1;
  row(ctr) = (i-1)*(Nx+2)+j;
  col(ctr) = (i-2)*(Nx+2)+j;
  A(ctr)   = 0.5*ds;
      
end
end

A = sparse(row,col,A);
