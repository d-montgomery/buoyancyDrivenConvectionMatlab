function [] = makePlots(U,V,T,t,xp,yp,XP,YP,XT,YT,XS,YS,ticks)

% Create Countour of PSI
PSI = streamFxn(U,V,xp,yp);
Svals = interp2(XP,YP,PSI,XS,YS);
Svals = Svals(:)';

hold off 

subplot(1,2,1)
contour(XP,YP,PSI, Svals,'-k'); 
xlabel('x [m]', 'FontSize',14)
ylabel('y [m]', 'FontSize',14)
title(['\psi(x,y,',num2str(t),')'],'FontSize',20)
xticks(ticks)
yticks(ticks)
axis([0 0.01 0 0.01])

subplot(1,2,2)
s = pcolor(XT,YT,T);
s.FaceColor = 'interp';
colorbar()
set(s, 'EdgeColor', 'none');
xlabel('x [m]')
ylabel('y [m]')
xticks(ticks)
yticks(ticks)
title('T_h')
title(['T(x,y,',num2str(t),')'],'FontSize',20)
hold on