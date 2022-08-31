function [ICx, ICy, ICp, ICt] = initial_conditions(Ly)

[~,~,~,~,Th,~,~,~] = parameters(Ly);

% Set Initial Conditions -------------------------------------------------
ICx = @(x,y) 0*x.*y;
ICy = @(x,y) 0*x.*y;
ICp = @(x,y) 0*x.*y;
ICt = @(x,y) Th + 0*x.*y ; % K
