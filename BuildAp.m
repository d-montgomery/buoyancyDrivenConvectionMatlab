function[Ap] = BuildAp(Nx,Ny,xp,yp,xu,yv)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% xp and yp is where p(x,y,t) is stored on grid
% xu is x-comp of where u(x,y,t) is stored on grid
% yv is y-comp of where v(x,y,t) is stored on grid
% mu is viscosity constant

% OUTPUT:
% A is a coefficient matrix

ctr = 0;

% Internal Nodes
for j = 1:Nx
for i = 1:Ny   

  xj = xp(j);
  yi = yp(i);
  dx = xu(j+1) - xu(j);
  dy = yv(i+1) - yv(i);
  
  if (j < Nx)
    xe = xp(j+1);
    ae = 1/(xe-xj)/dx;   
  else
    ae = 0;
  end
  
  if (j > 1)
    xw = xp(j-1);
    aw = 1/(xj-xw)/dx;
  else
    aw = 0;
  end
  
  if (i < Ny)
    yn = yp(i+1);
    an = 1/(yn-yi)/dy;
  else
    an = 0;
  end
  
  if (i > 1)
    ys = yp(i-1);
    as = 1/(yi-ys)/dy; 
  else 
    as = 0;
  end
  
  a =  -(ae + aw  + an + as);  
  
  ctr = ctr+1;
  row(ctr) = (i-1)*Nx + j;
  col(ctr) = (i-1)*Nx + j;
  A(ctr)   = a;
  
  if (j > 1)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-1)*Nx + j - 1;
    A(ctr)   = aw;
  end  
  
  if (j < Nx)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-1)*Nx + j + 1;
    A(ctr)   = ae;
  end
  
  if (i < Ny)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i)*Nx   + j;
    A(ctr)   = an;
  end
  
  if (i > 1)
    ctr = ctr+1;
    row(ctr) = (i-1)*Nx + j;
    col(ctr) = (i-2)*Nx + j;
    A(ctr)   = as;
  end  
  
  
end
end

% Add Lagrange Multiplier
for i = 1:Nx*Ny
  
  ctr = ctr+1;
  row(ctr) = i;
  col(ctr) = Nx*Ny + 1;
  A(ctr)   = 1.0;
  
end

ctr = ctr+1;
row(ctr) = Nx*Ny + 1;
col(ctr) = 1;
A(ctr)   = 1.0;

Ap = sparse(row,col,A);
