function [A] = BuildAu(Nx,Ny,xu,yu,xv,yv,mu,rho,dt)

% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% mu is viscosity constant
% rho is density
% dt is time step

% OUTPUT:
% A is a coefficient matrix

% Size of the A matrix
NNXY = (Nx+1)*(Ny+2);

% Initialize the counter used in building the sparse matrices
ctr = 0;

% Initialize row, col and A. S is the number of entries.
S = 4*(Nx - 2) + 4*(Ny + 2 - 1) + 5*(Nx -2) * (Ny+1 - 2);
row = zeros(1,S);
col = row;
A = row;

% Bottom BC
for j = 2:Nx
    % TS coefficient
    ctr = ctr + 1;
    row(ctr) = j;
    col(ctr) = j;
    A(ctr) = 1;
end

% Top BC 
for j = 2:Nx
    % TN coefficient
    ctr = ctr + 1;
    row(ctr) = (Ny+1)*(Nx+1) + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% Left and Right BC
for i = 1:Ny+2
    % Left BC
    % West
    ctr = ctr + 1;
    row(ctr) = (Nx+1)*(i-1)+1;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    % Right BC
    % East
    ctr = ctr + 1;
    row(ctr) = (Nx+1)*(i-1)+(Nx+1);
    col(ctr) = row(ctr);
    A(ctr) = 1;
end


% Fill Inner Entries of Matrices
for i = 2:Ny+1
    dy = yv(i)-yv(i-1);
    yn = yu(i+1);
    yp = yu(i);
    ys = yu(i-1);
    
    for j = 2:Nx
        dx = xv(j+1) - xv(j);
        xe = xu(j+1);
        xp = xu(j);
        xw = xu(j-1);
        
        % East
        ctr = ctr + 1;
        row(ctr) = (Nx+1)*(i-1)+j;
        col(ctr) = row(ctr) + 1;
        DE = - mu*dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = (Nx+1)*(i-1)+j;
        col(ctr) = row(ctr) + (Nx+1);
        DN = - mu*dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = (Nx+1)*(i-1)+j;
        col(ctr) = row(ctr) - 1;
        DW = - mu*dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = (Nx+1)*(i-1)+j;
        col(ctr) = row(ctr) - (Nx+1);
        DS = - mu*dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = (Nx+1)*(i-1)+j;
        col(ctr) = row(ctr);
        A(ctr) = -DE -DN -DW -DS + rho*dx*dy/dt;
    
    end
end
A = sparse(row, col, A, NNXY, NNXY);


