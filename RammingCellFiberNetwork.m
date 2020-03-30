% Simulation of a cell being pulled through the ECM matrix
clear;
close all;
addpath('movecellfunctions');
% Initialize the user-defined parameters (separate file)
AlexParameters;
% Perform a bunch of operations to initialize all variables
dt=dt0;
dthc=2*pi/Nc;
dthm=2*pi/Nm;
kc=kcsoft;
kbc = 0;
sc=[0:dthc:2*pi-dthc]';
sm=[0:dthm:2*pi-dthm]';
xm=rm*cos(sm);
ym=rm*sin(sm);
nucref=[xm ym];
xc=rc*cos(sc);
yc=rc*sin(sc);
x0=[xc yc];
x0_orig=x0;
hc=dthc*rc;
hm=dthm*rm;
bigR=1000;
% Make the ECM matrix
if (loadecm)
    load(ecmfile); % If there is a file to load
else
    [fiberlocs,connections]=makeECM(numfibers,rc,rfib,nodebound,eps,fibspace); % clean, new ECM
end
% The set up of the ECM is on the next 2 lines. The first function takes in
% the list of fiber centers and turns it into nFiber*nPtFiber vector of
% discretized fiber locations. ecmrefs is the intialized reference
% locations. The point of that is for the fibers to be initially unmoving.
% Because the lattice is random, every fiber gets a reference location that
% acts as an additional point the fiber is tethered to. In this way, the
% force is initially 0. 
ecmlocs=fiberlocs; 
ecmrefs=calcfibforce(ecmlocs,-1,connections,ecmlocs,rfib,0,0)+ecmlocs; 
pullingforce=zeros(Nc,2);
t=0;
tstar=t;
h=figure;
attached1=0;
attached2=0;
[prot1, prot2] = pickprots(Nc,ppoints(1),numfibers,Npfib,ecmlocs,xc,yc);
prot1=77;
prot2=9;
%prot1 = 61; % for sparse
%prot2 = 75; % for sparse
jPtslist=prot1;
nucpos=[0 0];
% Add the protrusions to the force vector pullingforce. 
pullingforce=pullingforce+addpinch(fmag,Nc,sc,prot1,xc,yc,rc)+addpinch(fmag,Nc,sc,prot2,xc,yc,rc);
retracting=0;
while (length(jPtslist) < length(ppoints))
    minnoded1=10;
    minnoded2=10;
    % attached means if the cell is attached or not. So if it's not
    % attached, check how close the end of the protrusion is to a fiber.
    % retracting is if the protrusion went too far and it's coming back in
    if (~attached1 && ~retracting)
        for ifibNum=1:numfibers
            for iN=(ifibNum-1)*Npfib+1:ifibNum*Npfib
                [d1,~] = min(sqrt((xc(prot1)-ecmlocs(iN,1)).^2+(yc(prot1)-ecmlocs(iN,2)).^2));
                if (d1 < minnoded1)
                    minnoded1=d1;
                    theNode1=iN;
                end
            end
        end
        tstar=t;
        % Check how long the protrusion is getting 
        protlen1=sqrt((xc(prot1)-mean(xc)).^2+(yc(prot1)-mean(yc)).^2)-rc;
    end 
    if (~attached2 && ~retracting)
        for ifibNum=1:numfibers
            for iN=(ifibNum-1)*Npfib+1:ifibNum*Npfib
                [d2,~] = min(sqrt((xc(prot2)-ecmlocs(iN,1)).^2+(yc(prot2)-ecmlocs(iN,2)).^2));
                if (d2 < minnoded2)
                    minnoded2=d2;
                    theNode2=iN;
                end
            end
        end
        tstar=t;
        % Check how long the protrusion is getting 
        protlen2=sqrt((xc(prot2)-mean(xc)).^2+(yc(prot2)-mean(yc)).^2)-rc;
    end 
    % Bring the protrusions back in if one of them gets too long.
    % No harm no foul, it just means there was no ECM fiber there.
    if (protlen1 > maxprotlen || protlen2 > maxprotlen) % retract a protrsion if it gets too long
        retracting=1;
        % move fiber
        if (attached1)
            % Move the cell eps back, then snap it back to where it was
            % before
            xc(prot1)=xc(prot1)+eps*fcel(prot1,1)/norm(fcel(prot1,:));
            yc(prot1)=yc(prot1)+eps*fcel(prot1,2)/norm(fcel(prot1,:));
            attached1=0;
        elseif (attached2)
            % Move the cell eps back, then snap it back to where it was
            % before
            xc(prot2)=xc(prot2)+eps*fcel(prot2,1)/norm(fcel(prot2,:));
            yc(prot2)=yc(prot2)+eps*fcel(prot2,2)/norm(fcel(prot2,:));
            attached2=0;
        end
        kc=kcrigid;
        kbc=kbcrigid;
        pullingforce=pullingforce*0;
    end
    if (retracting)
        dt=dt0/dtfactor;
        if (t-tstar > largerdttime)  % after sufficient time put the timestep back to what it was
            dt=dt0;
        end
    end
    % If the min node distance is less than 2*epsilon (1 epsilon for the
    % fiber and 1 for the cortex node), bind the cortex node.
    % Protrusion must be longer than some minimum length for this to
    % happen. 
    if ((minnoded1 < 2*eps && protlen1 > minprotlen) || attached1)
        if (attached1==0)
            pullingforce=pullingforce-addpinch(fmag,Nc,sc,prot1,xc,yc,rc);
        end
        attached1=1;
        xc(prot1)=ecmlocs(theNode1,1); % move the cortex points onto the fiber
        yc(prot1)=ecmlocs(theNode1,2);
    end
    if ((minnoded2 < 2*eps && protlen2 > minprotlen) || attached2)
        if (attached2==0)
            pullingforce=pullingforce-addpinch(fmag,Nc,sc,prot2,xc,yc,rc);
        end
        attached2=1;
        xc(prot2)=ecmlocs(theNode2,1); % move the cortex points onto the fiber
        yc(prot2)=ecmlocs(theNode2,2);
    end
    if (attached1 && attached2)
        dt=dt0/dtfactor;             % shorten up the timestep because the cell is stiffer
        if (t-tstar > largerdttime)  % after sufficient time put the timestep back to what it was
            dt=dt0;
        end
    end
    % Routine force calculations - nucleus and cortex elastic force and
    % total force. 
    if (attached1 && attached2)
        % You have to remesh every once in a while! Here I'm doing it every
        % 0.1
        if (mod(t-tstar,0.01) < dt && t-tstar > dt)
            [newpos, x0,prot1,prot2,~,~] = remeshCortex(Nc,[xc yc],x0,prot1,prot2);
            [newnuc,nucref]=remeshNucleus(Nm,[xm ym],nucref);
            xc=newpos(:,1);
            yc=newpos(:,2);
            xm=newnuc(:,1);
            ym=newnuc(:,2);
        end
        fcel=calcnonunifelf(rc,xc,yc,x0(:,1),x0(:,2),gamc,gamc,kcsoft,kcrigid,prot1,prot2);
        fcel = fcel+calcnonunifbend(xc,yc,rc,x0,kbcrigid,prot1,prot2);
    else
        fcel=calcelasticforce(rc,xc,yc,gamc,kc);
        fcel = fcel+calcbendingforce(xc,yc,rc,kbc);
    end
    fmel=calcelasticforcerefs(rm,xm,ym,gamm,km,nucref);
    fmbend=calcbendingforcerefs(xm,ym,rm,kbm,nucref);
    fmemtot=(fmel+fmbend)*hm; %Force/length
    fcortot=(fcel+pullingforce)*hc;
    if (rigidecm) % If the ECM is rigid, separate set of steps. Enforce the 0 velocity boundary 
                  % condition on the fibers, then solve the system for the
                  % motion of the cell.
        alllocs=[xm ym; xc yc; ecmlocs];
        M=formMatrix(alllocs,eps,mu);
        f=zeros(Nc+Nm,1);
        for ik=1:Nm
            f(2*ik-1)=fmemtot(ik,1);
            f(2*ik)=fmemtot(ik,2);
        end
        for ik=1:Nc
            f(2*ik-1+2*Nm)=fcortot(ik,1);
            f(2*ik+2*Nm)=fcortot(ik,2);
        end
        wouldbeus=M(2*(Nc+Nm)+1:length(M),1:2*(Nc+Nm))*f;
        fibforce1 = M(2*(Nc+Nm)+1:length(M),2*(Nc+Nm)+1:length(M))\(-wouldbeus);
        totforce1=[f; fibforce1];
        uvalsv=M*totforce1;
        for ik=1:length(uvalsv)/2
            uvals(ik,1)=uvalsv(2*ik-1);
            uvals(ik,2)=uvalsv(2*ik);
        end
        for ik=1:length(fibforce1)/2
            fibforce(ik,1)=fibforce1(2*ik-1);
            fibforce(ik,2)=fibforce1(2*ik);
        end
    else % otherwise it's what you would expect. Calculate the fiber force, append to the list, solve
        fibforce=calcfibforce(ecmlocs,kteth,connections,ecmrefs,rfib,gamf,kf);
        alllocs=[xm ym; xc yc; ecmlocs];
        totforce=[fmemtot; fcortot; fibforce];
        M=formMatrix(alllocs,eps,mu);
        f=zeros(length(alllocs)*2,1);
        for ik=1:length(alllocs)
            f(2*ik-1)=totforce(ik,1);
            f(2*ik)=totforce(ik,2);
        end 
        uvals1=M*f;
        for ik=1:length(alllocs)
            uvals(ik,1)=uvals1(2*ik-1);
            uvals(ik,2)=uvals1(2*ik);
        end 
        uvals=uvals+(log(bigR)-0.5)/(4*pi*mu)*sum(totforce); % add the constraint velocity from 
        % our paper (average velocity 0 on the "large circle")
    end
    if (mod(floor(t/dt+1e-5),saveEvery)==0)
        plotmigcell; % Make a plot of what's happening
    end
    % Advance the locations
    alllocs=alllocs+dt*uvals;
    xm=alllocs(1:Nm,1);
    ym=alllocs(1:Nm,2);
    xc=alllocs(Nm+1:Nm+Nc,1);
    yc=alllocs(Nm+1:Nm+Nc,2);
    ecmlocs=alllocs(Nm+Nc+1:end,:);
    t=t+dt;
    % If the whole system has completely settled down, restart. 
    if (max(sqrt(uvals(:,1).^2+uvals(:,2).^2)) < eps && retracting)
        retracting=0;
        dt=dt0;
        tstar=t; % start over
        kc=kcsoft; % back to a soft cortex
        kbc = 0;
        nucpos=[nucpos; mean(xm) mean(ym)];
        [prot1, prot2] = pickprots(Nc,ppoints(length(jPtslist)+1),numfibers,Npfib,ecmlocs,xc,yc);
        prot1=71;
        prot2=7;
        jPtslist=[jPtslist prot1];
        pullingforce=pullingforce+addpinch(fmag,Nc,sc,prot1,xc,yc,rc)+...
            addpinch(fmag,Nc,sc,prot2,xc,yc,rc);   % add the pinch and start over
    end
    % Take off the attachments if the cell has settled down. 
    if (max(sqrt(uvals(:,1).^2+uvals(:,2).^2)) < eps && attached1 && attached2)
        % move fiber
        thefib1=ceil(theNode1/Npfib);
        thefib2=ceil(theNode2/Npfib);
        [newpos, ~,prot1,prot2,np1,np2] = remeshCortex(Nc,[xc yc],x0,prot1,prot2);
        xc=newpos(:,1);
        yc=newpos(:,2);
        % Move the ECM fiber 2*eps away (like the cell "spitting it out")
        ecmlocs((thefib1-1)*Npfib+1:thefib1*Npfib,:)=ecmlocs((thefib1-1)*Npfib+1:thefib1*Npfib,:)+...
            2*eps*np1;
        ecmlocs((thefib2-1)*Npfib+1:thefib2*Npfib,:)=ecmlocs((thefib2-1)*Npfib+1:thefib2*Npfib,:)+...
            2*eps*np2;
        attached1=0;
        attached2=0;
        retracting=1;
        x0=x0_orig;
    end
    if(prot1==prot2)
        keyboard
    end
end
% After completion, determine the total amount the nucleus moved
normd=norm(nucpos(end,:))
protds=nucpos-[0 0; nucpos(1:length(ppoints)-1,:)];
protds(:,3)=sqrt(protds(:,1).^2+protds(:,2).^2);
% Number of nucleus nodes ahead of ecm 32 and 39
xline = (ecmlocs(39,1)-ecmlocs(32,1))/(ecmlocs(39,2)-ecmlocs(32,2))*(ym-ecmlocs(32,2))+ecmlocs(32,1);
nucpen=sum(xm > xline)/Nm
corpen=sum(xc > xline)/Nc



