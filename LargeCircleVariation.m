% Test for assumption of constant velocity
set(0,'DefaultAxesFontName','Times New Roman')
N=1000;
h=2*pi/N;
s=[0:h:2*pi-h]';
eps=0;
mu=1;
for logbigR=1:20
    bigR=10^(logbigR);
    xl=bigR*cos(s);
    yl=bigR*sin(s);
    f=[1 0];
    [~,wouldbeus]=regStokesSolve(xl,yl,0,0,f,eps,mu);
    v_r(log10(bigR))=(max(wouldbeus(:,1))-mean(wouldbeus(:,1)))/mean(wouldbeus(:,1));
end 
plot(1:20,abs(v_r)*100,'bo')
set(gca,'FontSize',16)
xlabel('$R$','Interpreter','latex')
ylabel({'$\sigma_u(R)$ (\%)'},'Interpreter','latex')
xticks([1:20])
xticklabels({'10^{1}','','','','10^{5}','','','','','10^{10}','','','','','10^{15}','','','','','10^{20}'})
yticks([0:5:60])
yticklabels({'','','10','','20','','30','','40','','50','','60'})
xlim([0 21])
hold on
plot(0:21,10*ones(1,22),':k')
plot(0:21,5*ones(1,22),':k')



