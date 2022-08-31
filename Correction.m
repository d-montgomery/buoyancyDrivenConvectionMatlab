function[U,V,P] = Correction(U,V,P,PHI,Nx,Ny,rho,dt,xp,yp)
% INPUTS:
% U, V, P, are the velocities and pressure
% PHI is the pseudo pressure
% Nx and Ny are the number of FV cells in each direction
% rho is density
% dt is time step
% xp and yp is where p(x,y) is stored on grid


% OUTPUT: Corrected U,V,P

% update Pressure
P = P + PHI;

% update U
for i = 2:Ny+1
for j = 2:Nx

U(i,j) = U(i,j) - dt/rho*( PHI(i-1,j) - PHI(i-1,j-1) )/( xp(j)-xp(j-1) );

end
end

% update V
for i = 2:Ny
for j = 2:Nx+1

V(i,j) = V(i,j) - dt/rho*( PHI(i,j-1) - PHI(i-1,j-1) )/( yp(i)-yp(i-1) );

end
end





