% Make the pressure/velocity plots in Fig. 5
xmax=1.6;
ymin=-0.5;
x=[-0.5:0.02:xmax];
y=[ymin:0.02:0.6];
[targsX,targsY]=meshgrid(x,y);
[p,u]=regStokesSolve(targsX(:),targsY(:),alllocs(:,1),alllocs(:,2),totforce,eps,mu);
uextra=(log(bigR)-0.5)/(4*pi*mu)*sum(totforce);
p=reshape(p,length(y),length(x));
ux=reshape(u(:,1),length(y),length(x))+uextra(1);
uy=reshape(u(:,2),length(y),length(x))+uextra(2);
figure;
imagesc([-0.5 xmax],[ymin 0.6],p)
set(gca,'YDir','Normal')
colormap(flipud(hot))
hold on
plotmigcell;
xlim([-0.5 xmax])
ylim([ymin 0.6])
pbaspect([0.5+xmax 0.6-ymin 1])
colorbar
figure;
imagesc([-0.5 xmax],[ymin 0.6],ux)
set(gca,'YDir','Normal')
colormap(flipud(hot))
colorbar
hold on
plotmigcell;
xlim([-0.5 xmax])
ylim([ymin 0.6])
pbaspect([0.5+xmax 0.6-ymin 1])
colorbar
figure;
ev=5;
plotmigcell;
hold on
quiver(targsX(1:ev:end,1:ev:end),targsY(1:ev:end,1:ev:end),ux(1:ev:end,1:ev:end),uy(1:ev:end,1:ev:end),...
    'k','LineWidth',1.5)
xlim([-0.5 xmax])
ylim([ymin 0.6])
pbaspect([0.5+xmax 0.6-ymin 1])