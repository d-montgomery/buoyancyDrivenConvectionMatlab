function [xp,yp,xt,yt,xu,yu,xv,yv] = BuildGrid(Nx,Ny,Lx,Ly,GRID,PLOT)
% INPUTS:
% Nx and Ny are the number of FV cells in each direction
% Lx and Ly are the length of the domain in the x/y directions
% GRID = 0 or 1 to specify uniform vs. Chebshev grid
% PLOT = 0 or 1 to specify if a plot of the grid is desired

% OUTPUTS:
% Returns vectors of where each variable is stored on FV grid
% e.g. p(x_j, y_i) is located at xp(j), yp(i)

if GRID == 1  % Chebshev Grid
    x = -Lx / 2*cos((0:Nx)*pi/Nx) +  Lx / 2;
    y = -Ly / 2*cos((0:Ny)*pi/Ny) +  Ly / 2;
else % Uniform Grid
    dx = Lx / Nx;
    dy = Ly / Ny;
    x = 0:dx:Lxb;
    y = 0:dy:Lyb;
end

% Length of each cell (for nonuniform grid)
dx = x(2:Nx+1) - x(1:Nx);
dy = y(2:Ny+1) - y(1:Ny);

% Location of pressure (mid points of cells)
xp = x(1:Nx) + dx/2;
yp = y(1:Ny) + dy/2;

% Location of u
xu = x;
yu = [y(1), y(1:Ny) + dy/2, y(Ny+1)];

% Location of v
xv = [x(1), x(1:Nx) + dx/2, x(Nx+1)];
yv = y; 

% Location of Temperature
xt = xv;
yt = yu;

if PLOT == 1
    % Marker size
    M = 15;
    
    [Px,Py] = meshgrid(xp,yp);
    [Ux,Uy] = meshgrid(xu,yu);
    [Vx,Vy] = meshgrid(xv,yv);
    [Tx,Ty] = meshgrid(xt,yt);
    
    fig1 = figure(1);
    pos_fig1 = [0 0 750 600];
    set(fig1,'Position',pos_fig1)
    hold on
    axis([x(1), x(Nx+1), y(1), y(Ny+1)])
    plot(Px, Py,'k.','MarkerSize', M)
    plot(Ux, Uy,'bs','MarkerSize', M)
    plot(Vx, Vy, 'b^','MarkerSize', M)
    plot(Tx, Ty, 'ro','MarkerSize', M)
    
    % Plot vertical and horizontal grid lines
    for i = 1: Nx+1
        plot([x(i), x(i)], [y(1), y(Ny+1)], 'r')
    end
    for i = 1:Ny+1
        plot([x(1), x(Nx+1)], [y(i), y(i)], 'r')
    end
    
    % Labels, legend, etc
    fs = 16;
    xlabel('x', 'FontSize', fs)
    ylabel('y', 'FontSize', fs)
      
    
    % for legend puproses
    l1 = plot(nan, nan,'k.','MarkerSize', M, 'DisplayName', 'p');
    l2 = plot(nan, nan,'bs', 'MarkerSize', M, 'DisplayName', 'u');
    l3 = plot(nan, nan, 'b^', 'MarkerSize', M, 'DisplayName', 'v');
    l4 = plot(nan, nan, 'ro', 'MarkerSize', M, 'DisplayName', 'T');
    legend([l1,l2,l3,l4], 'FontSize', fs,'Location','northeastoutside')
    
end
end