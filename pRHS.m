function[F] = pRHS(Nx,Ny,U,V,rho,dt,xu,yv)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% U is predictive velocity U
% V is predictive velocity V
% rho is density
% dt is time step
% xu is x-comp of where u(x,y) is stored on grid
% yv is y-comp of where v(x,y) is stored on grid

% OUTPUT:
% F is a vector

F = zeros(Ny*Nx+1,1);

for i = 1:Ny
for j = 1:Nx
  
    dx = xu(j+1) - xu(j);  
    dy = yv(i+1) - yv(i);
   
    F( (i-1)*Nx + j ) = rho/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx +  ( V(i+1,j+1) - V(i,j+1) )/dy );
    
end  
end

F(Ny*Nx+1) = 0;





