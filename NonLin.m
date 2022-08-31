function[NLU,NLV,NLT] = NonLin(U,V,T,Nx,Ny,rho,xu,yu,xv,yv)
% INPUTS:
% U, V, T, are the velocities and temp
% Nx and Ny are the number of FV cells in each direction
% rho is density
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where u(x,y) is stored on grid

% OUTPUT: Matrices containing the nonlinear terms for U,V,T

% --- Compute NLU for x-momentum equation --------------------------------

NLU = zeros(Ny+1,Nx);

for i = 2:Ny+1
for j = 2:Nx
    dy = yv(i) - yv(i-1);

    % Top and Bottom Edges of cell
    ys = yv(i-1); yn = yv(i);

    % Location of Nodes
    XE = xu(j+1);  XP = xu(j);  XW = xu(j-1);
    YN = yu(i+1);  YP = yu(i);  YS = yu(i-1);


    % East
    ue = 0.5*(U(i,j+1) + U(i,j));
    me = 0.5 * rho * dy * (U(i,j) + U(i,j+1));

    % North
    if i == Ny+1 
        un = U(i+1,j); % Don't interpolate at boundary
    else
        un = U(i,j) + (yn - YP) * ( U(i+1,j) - U(i,j) ) / ( YN - YP );
    end
    mn = 0.5 * rho * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


    % West
    uw = 0.5*(U(i,j) + U(i,j-1));
    mw = 0.5 * rho *dy * (U(i,j) + U(i,j-1));

    % South
    if i == 2
        us = U(i-1,j); % Don't interpolate at boundary
    else
        us = U(i,j) + (ys-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
    end
    ms = 0.5 * rho * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

    NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
end
end

% --- Compute NLV for y-momentum equation --------------------------------
NLV = zeros(Ny,Nx+1);


for i = 2:Ny
for j = 2:Nx+1

    % Cell height/width
    dx = xu(j) - xu(j-1);

    % Left and Right Edges of cell
    xw = xu(j-1); xe = xu(j);

    % Location of Nodes
    XE = xv(j+1);  XP = xv(j);  XW = xv(j-1);
    YN = yv(i+1);  YP = yv(i);  YS = yv(i-1);

    % East
    if j == Nx+1
        ve = V(i,j+1); % don't interpolate on boundary
    else
        ve = V(i,j) + (xe-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
    end
    me = 0.5*rho*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

    % North
    vn = 0.5*( V(i+1,j) + V(i,j) );
    mn = 0.5*rho* dx * ( V(i+1,j) + V(i,j) );


    % West
    if j == 2
        vw = V(i,1); % don't interpolate on boundary
    else
        vw = V(i,j) + (xw-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
    end
    mw = 0.5*rho*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

    % South
    vs = 0.5*( V(i,j) + V(i-1,j) );
    ms = 0.5*rho*dx * ( V(i,j) + V(i-1,j) );

    NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
end
end

% --- Compute NLT for Heat Equation --------------------------------
NLT=zeros(Ny+2,Nx+2);

% Internal Nodes
for j = 2:Nx+1
for i = 2:Ny+1   

  % compute size of cell faces
  dx = xu(j) - xu(j-1);  
  dy = yv(i) - yv(i-1);
  
  % compute fluid velocites on four faces
  Uw = U(i,j-1);
  Ue = U(i,j);
  Vs = V(i-1,j);
  Vn = V(i,j);
  
  % east face
  if (j == Nx+1)
    ae = Ue*dy*T(i,j+1);
  else
    ae = 0.5*Ue*dy*( T(i,j) + T(i,j+1) );
  end
  
  if (j == 2)
    aw = -Uw*dy*T(i,j-1);
  else
    aw = -0.5*Uw*dy*( T(i,j) + T(i,j-1) );
  end
  
  if (i == 2)
    as = -Vs*dx*T(i-1,j);
  else
    as = -0.5*Vs*dx*( T(i,j)+T(i-1,j) );
  end
  
  if (i == Ny+1)
    an = Vn*dx*T(i+1,j);
  else
    an = 0.5*Vn*dx*( T(i,j)+T(i+1,j) );
  end
  
  NLT(i,j) = ae + an + aw + as;  
      
end
end