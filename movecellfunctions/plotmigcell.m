% Script to make a movie of the cell moving
if (t<dt)
    quiver(xc,yc,fcortot(:,1),fcortot(:,2),2,'k','LineWidth',1.25)
    hold on
end
%quiver(xm,ym,fmemtot(:,1),fmemtot(:,2))
%quiver(alllocs(2*N+1:end,1),alllocs(2*N+1:end,2),totforce(2*N+1:end,1),totforce(2*N+1:end,2),...
%    'AutoScaleFactor',0.5)
plot([xm; xm(1)],[ym; ym(1)],'m','LineWidth',2)
hold off
hold on
plot([xc; xc(1)],[yc; yc(1)],'color',[0 0.7 0.3],'LineWidth',2)
%scatter(alllocs(Nm+Nc+1:end,1),alllocs(Nm+Nc+1:end,2),72,'b','filled')
scatter(alllocs(Nm+Nc+1:end,1),alllocs(Nm+Nc+1:end,2),'b','filled')
%scatter(fiberlocs([32 39],1),fiberlocs([32 39],2),108,'kx','LineWidth',2.5)
xlim([-rc-0.1 nodebound+0.1])
ylim([-nodebound/2-0.1 nodebound/2+0.1])
% xlim([-0.6 2])
% ylim([-1 1])
% pbaspect([2.6 2 1])
pbaspect([nodebound+rc nodebound 1])
xlabel({'$x$'},'Interpreter','latex')    
ylabel({'$y$'},'Interpreter','latex')
str=sprintf('$t=$ %1.2f',t);
title(str,'Interpreter','latex')
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',16)
hold off
if (t==0)
    movieframes=getframe(h);
elseif (mod(floor(t/dt+1e-5),saveEvery)==0)
    movieframes(length(movieframes)+1)=getframe(h);
end