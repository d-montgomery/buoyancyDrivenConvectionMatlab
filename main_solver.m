function [T] = main_solver(Nx,Ny,Lx,Ly,Tf,dt,flags)
% Inputs:
% Nx and Ny are the number of FV cells in x/y directions 
% Lx and Ly are the length of the domain in the x/y directions
% Tf is the final time
% dt is the temporal step size
% flags is an array with flags = {plotFlag, saveFlag, gridFlag, incr};

% unpack flags
plotFlag = flags(1);
saveFlag = flags(2);
gridFlag = flags(3);
incr = flags(4);

% Get all necessary physical parameters
[mu,rho,alpha,~,~,~,~,~] = parameters(Ly);

% Get grid for FVM 
[xp,yp,xt,yt,xu,yu,xv,yv] = BuildGrid(Nx,Ny,Lx,Ly,gridFlag,plotFlag);
[XP, YP] = meshgrid(xp,yp);
[XT, YT] = meshgrid(xt,yt);
[XU, YU] = meshgrid(xu,yu);
[XV, YV] = meshgrid(xv,yv);

% Calculate number of time steps Nt
t = 0 : dt : Tf;
Nt = length(t) - 1;

% Unpack Boundary Condition Info
[at,bt,BC] = boundary_conditions(Ly);

% Get Matrices Au, Av, and At
Au = BuildAu(Nx,Ny,xu,yu,xv,yv,mu,rho,dt);
Av = BuildAv(Nx,Ny,xu,yu,xv,yv,mu,rho,dt);
At = BuildAt(Nx,Ny,xt,yt,xu,yv,alpha,dt,at,bt);

% Get Pressure Matrix and initialize PHI
Ap = BuildAp(Nx,Ny,xp,yp,xu,yv);
PHI = zeros(Nx,Ny);

% Store factorizations for efficient solving
Au = decomposition(Au);
Av = decomposition(Av);
Ap = decomposition(Ap);
At = decomposition(At);

% Initial Conditions
[ICx, ICy, ICp, ICt] = initial_conditions(Ly);
U0 = ICx(XU,YU);
V0 = ICy(XV,YV);
P0 = ICp(XP,YP);
T0 = ICt(XT,YT);

% Forcing Terms at t = t0
T0v = interp_T_to_v(Ny,yu,T0); % Interp T0 to v momentum cells
[Qxo,Qyo,Qto] = body_forcing(XU,XT,T0v,Ly);

% --- One step of Forward Euler ---
[NLU,NLV,NLT] = NonLin(U0,V0,T0,Nx,Ny,rho,xu,yu,xv,yv);
NLU = 2/3*NLU;
NLV = 2/3*NLV;
NLT = 2/3*NLT;
NLU0 = 0*NLU;
NLV0 = 0*NLV;
NLT0 = 0*NLT;


% Initiate Stream Line plot
if plotFlag == 2
    % Components for Stream Fxn Figure
    ticks = 0:0.002:0.01;
    NSX=10;
    NSY=10;
    ys = linspace(min(yp),max(yp),NSX);
    xs = linspace(min(xp),max(xp),NSY);
    [XS,YS] = meshgrid(xs,ys);
    figure(1)
    if (saveFlag==1)
        set(gcf,'position',[30,30,1400,800])
        vid = VideoWriter(['vidTf',num2str(Tf),'.mp4'],'MPEG-4');
        open(vid);
    end
end

% Step in time
for k = 2:Nt+1
    % --- Solve for Temp ---
    % Body Forcing terms for temp
    [~,~,Qt] = body_forcing(XU,XT,T0v,Ly);
    
    % Boundary Conditions
    T1 = BC{3,1}(yt',t(k)); T2 = BC{3,2}(yt',t(k));
    T3 = BC{3,3}(xt,t(k));  T4 = BC{3,4}(xt,t(k));
    BCt = {T1,T2,T3,T4};
    
    Ft = tRHS(Nx,Ny,T0,NLT,NLT0,Qt,Qto,xt,yt,xu,yv,dt,alpha,BCt);
    T = At \ Ft;
    T = reshape(T, Nx+2,Ny+2)';
    
    % --- Solve for U,V,P ----
    % Body forcing 
    Tv = interp_T_to_v(Ny,yu,T); % Interp T to v-momentum cells
    [Qx,Qy,~] = body_forcing(XU,XT,Tv,Ly);
    
    % Boundary Conditions for u and v 
    g1 = BC{1,1}(yu',t(k)); g2 = BC{1,2}(yu',t(k)); 
    g3 = BC{1,3}(xu,t(k));  g4 = BC{1,4}(xu,t(k)); 
    h1 = BC{2,1}(yv',t(k)); h2 = BC{2,2}(yv',t(k)); 
    h3 = BC{2,3}(xv,t(k));  h4 = BC{2,4}(xv,t(k));
    
    % Pass BC's in Cell Array
    BCx = {g1,g2,g3,g4};    BCy = {h1,h2,h3,h4}; 
    
    % Prediction - Solve x-equation
    Fu = uRHS(Nx,Ny,U0,NLU,NLU0,P0,Qx,Qxo,dt,rho,mu,xu,yu,xv,yv,BCx);
    U = Au \ Fu; 
    U = reshape(U, Nx+1,Ny+2)';

    % Prediction - Solve y-equation
    Fv = vRHS(Nx,Ny,V0,NLV,NLV0,P0,Qy,Qyo,dt,rho,mu,xu,yu,xv,yv,BCy);
    V = Av \ Fv;
    V = reshape(V, Nx+2,Ny+1)';

    % Correction - Solve PHI Equation
    Fp = pRHS(Nx,Ny,U,V,rho,dt,xu,yv);
    TMPP = Ap\Fp;
    for i=1:Ny
        PHI(i,:) = TMPP( (i-1)*(Nx) + 1 : i*(Nx));
    end 
    
    % Correction - Get Current U, V and P
    [U,V,P] = Correction(U,V,P0,PHI,Nx,Ny,rho,dt,xp,yp);
    
    % Update Velocities
    U0 = U;
    V0 = V;
    P0 = P;
    T0 = T;
    Qxo = Qx;
    Qyo = Qy;
    Qto = Qt;
    NLU0 = NLU;
    NLV0 = NLV;
    NLT0 = NLT;
    [NLU,NLV,NLT] = NonLin(U,V,T,Nx,Ny,rho,xu,yu,xv,yv);
    

    % plot solution at various values of time
    if plotFlag == 2 && ~mod(k,incr)
        
        makePlots(U,V,T,t(k),xp,yp,XP,YP,XT,YT,XS,YS,ticks)
        pause(0.1)
        
        if (saveFlag==1) % Save frame for video
            frame = getframe(gcf);
            writeVideo(vid,frame);
        end
    end
end

% Save a movie of Stream Lines
if (plotFlag == 2) && (saveFlag==1)
    % Output the movie as an avi file
    close(vid);
end