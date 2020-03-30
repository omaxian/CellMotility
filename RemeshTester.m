close all;
addpath('movecellfunctions')
Nc=25;
rc=1;
prot1=floor(Nc*0.75);
prot2=ceil(Nc*0.25);
h=2*pi*rc/Nc;
sc=[0:h/rc:2*pi-h/rc]';
xc=rc*cos(sc);
yc=rc*sin(sc);
ksoft=1;
khard=100;
gamsoft=0;
gamhard=0;
t=0;
dt=0.01;
t_final=8;
pos=[xc yc];
x_new_ref=pos;
g=figure;
uvals1=1;
disp=[pos(prot1,1)];
while (t < t_final+dt)
    if (t > t_final-1-dt && t < t_final-1+1e-10)
        posb4=pos;
        [pos, x_new_ref,prot1,prot2,np1,np2] = remeshCortex(Nc,pos,x_new_ref,prot1,prot2);
    end
    [force,pts4ten,ten]=meshCheckerForce(rc,pos(:,1),pos(:,2),x_new_ref(:,1),x_new_ref(:,2),...
        gamsoft,gamhard,ksoft,khard,prot1,prot2);
    force = force+h*calcbendingforcerefs(pos(:,1),pos(:,2),rc,0.1,x_new_ref);
%     sum(force)
    M=formMatrix(pos,h,1);
    f=zeros(length(pos)*2,1);
    for ik=1:length(pos)
        f(2*ik-1)=force(ik,1);
        f(2*ik)=force(ik,2);
    end 
    uvals1=M*f;
    uvals=zeros(Nc,2);
    for ik=1:length(pos)
        uvals(ik,1)=uvals1(2*ik-1);
        uvals(ik,2)=uvals1(2*ik);
    end 
%     uvals(Nc/2,1)
    plot([pos(prot1:end,1); pos(1:prot2,1)],[pos(prot1:end,2); pos(1:prot2,2)],'-bo','LineWidth',0.85);
    hold on
    plot(pos(prot2:prot1,1),pos(prot2:prot1,2),'-rd','LineWidth',0.85);
%     scatter(pts4ten(:,1),pts4ten(:,2),60,ten,'filled')
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
%     legend({'$k=1$ pN/$\mu$m, $r=1$ $\mu$m','$k=100$ pN/$\mu$m, $r=0.5$ $\mu$m'},'Interpreter','latex')
%     colorbar;
%     caxis([0 60])
    xlim([-1 2])
    ylim([-1.5 1.5])
    pbaspect([1 1 1])
    hold off
    movieframes(floor(t/dt+1e-5)+1)=getframe(g);
    pos=pos+dt*uvals;
    t=t+dt;
    disp=[disp pos(prot1,1)];
end

function addforce=calcp0(xc,yc,force,rc)
    % Calculate normals
    normals=zeros(length(xc),2);
    tsum=0;
    N=length(xc);
    h=2*pi*rc/N;
    for iPt=1:N
        indexm1=iPt-1;
        if (indexm1==0)
            indexm1=N;
        end 
        indexp1=iPt+1;
        if (indexp1==N+1)
            indexp1=1;
        end 
        tao=[xc(indexp1)-xc(indexm1); yc(indexp1)-yc(indexm1)]/(2*h);
        tao=tao/norm(tao);
        % Rotate 90 degrees
        normals(iPt,:)=([0 -1; 1 0]*tao)';
        tsum=tsum-dot(force(iPt,:),normals(iPt,:));
    end
    tsum=tsum*h;
    p0=tsum/(2*pi*rc);
    addforce=normals*p0;
end
    
    