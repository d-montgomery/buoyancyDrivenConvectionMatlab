function PSI = streamFxn(U,V,xp,yp)

Nx = length(xp);
Ny = length(yp);
PSI = zeros(Ny, Nx);

for i = 2:Ny
    dy = yp(i) - yp(i-1);
    un = 0.5 * ( U(i+1,1) + U(i+1,2) );
    up = 0.5 * ( U(i,1) + U(i,2) );
    PSI(i,1) = PSI(i-1,1) +  0.5 * dy *( un + up );
end

for j = 2:Nx
    dx = xp(j) - xp(j-1);
    ve = 0.5 * ( V(2:Ny+1 , j+1) + V(1:Ny , j+1) );
    vp = 0.5 * ( V(2:Ny+1 , j) + V(1:Ny , j) );
    PSI(1:Ny,j) = PSI(1:Ny,j-1) - 0.5 * dx .* ( ve + vp );
end

    
    
    
    
    
    
