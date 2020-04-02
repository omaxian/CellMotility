% Simulation of a cell being rammed through the ECM matrix
% This corresponds to mechanism 2 in the paper. 
addpath('movecellfunctions');
close all;
% Initialize the user-defined parameters (separate file)
CellParameters_Mech2;
% Perform a bunch of operations to initialize all variables
dt=dt0;
dthc=2*pi/Nc;
dthm=2*pi/Nm;
kc=kcsoft;
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
t0=1;
% Make the ECM matrix
if (loadecm)
    load(ecmfile); % If there is a file to load
else
    [fiberlocs,connections]=makeECM(numfibers,rc,nodebound,eps,fibspace); % clean, new ECM
end
% The set up of the ECM is on the next 2 lines. ecmrefs is the intialized reference
% locations. The point of that is for the fibers to be initially unmoving.
% Because the lattice is random, every fiber gets a reference location that
% acts as an additional point the fiber is tethered to. In this way, the
% force is initially 0. 
ecmlocs=fiberlocs; 
ecmrefs=calcfibforce(ecmlocs,-1,connections,ecmlocs)+ecmlocs; 
pullingforce=zeros(Nc,2);
t=0;
tstar=t;
h=figure;
attached1=0;
attached2=0;
[prot1, prot2] = pickprots(Nc,ppoints(:,1),t0,km);
jPtslist=[prot1; prot2];
nucpos=[0 0];
% Add the protrusions to the force vector pullingforce. 
pullingforce=pullingforce+addpinch(fmag,Nc,sc,prot1)+addpinch(fmag,Nc,sc,prot2);
retracting=0;
AcorError = 0;
AnucError = 0;
CorAR = 1;
theNode2 = [];
while (length(jPtslist(1,:)) < length(ppoints(1,:)))
    minnoded1=10;
    minnoded2=10;
    % attached means if the cell is attached or not. So if it's not
    % attached, check how close the end of the protrusion is to a fiber.
    % retracting is if the protrusion went too far and it's coming back in
    if (~attached1 && ~retracting)
        dists = sqrt(sum(([xc(prot1) yc(prot1)]-ecmlocs) .*([xc(prot1) yc(prot1)]-ecmlocs),2));
        dists(theNode2) = inf; % prevents nodes from binding to the same site
        [minnoded1,theNode1] = min(dists);
        tstar=t;
        % Check how long the protrusion is getting 
        protlen1=sqrt((xc(prot1)-mean(xc)).^2+(yc(prot1)-mean(yc)).^2)-rc;
    end 
    if (~attached2 && ~retracting)
        dists = sqrt(sum(([xc(prot2) yc(prot2)]-ecmlocs) .*([xc(prot2) yc(prot2)]-ecmlocs),2));
        dists(theNode1) = inf; % prevents nodes from binding to the same site
        [minnoded2,theNode2] = min(dists);
        tstar=t;
        % Check how long the protrusion is getting 
        protlen2=sqrt((xc(prot2)-mean(xc)).^2+(yc(prot2)-mean(yc)).^2)-rc;
    end 
    % Bring the protrusions back in if one of them gets too long (no ECM contact)
    if (~retracting && ~(attached1 && attached2) && ...
        (protlen1 > maxprotlen || protlen2 > maxprotlen))
        retracting=1;
        if (attached1)
            xc(prot1) = xc(prot1)+eps*fcel(prot1,1)/norm(fcel(prot1,:));
            yc(prot1) = yc(prot1)+eps*fcel(prot1,2)/norm(fcel(prot1,:));
        end
        if (attached2)
            xc(prot2) = xc(prot2)+eps*fcel(prot2,1)/norm(fcel(prot2,:));
            yc(prot2) = yc(prot2)+eps*fcel(prot2,2)/norm(fcel(prot2,:));
        end
        attached1 = 0;
        attached2 = 0;
        kc=kcrigid;
        kbc=kbcrigid;
        pullingforce=pullingforce*0;
    end
    % If the min node distance is less than 2*epsilon (1 epsilon for the
    % fiber and 1 for the cortex node), bind the cortex node.
    % Protrusion must be longer than some minimum length for this to
    % happen. 
    if (minnoded1 < 2*eps || attached1)
        if (attached1==0) % take off the pulling force associated with prot1
            pullingforce=pullingforce-addpinch(fmag,Nc,sc,prot1);
        end
        attached1=1;
        xc(prot1)=ecmlocs(theNode1,1); % move the cortex points onto the fiber
        yc(prot1)=ecmlocs(theNode1,2);
    end
    if (minnoded2 < 2*eps || attached2)
        if (attached2==0) % take off the pulling force associated with prot2
            pullingforce=pullingforce-addpinch(fmag,Nc,sc,prot2);
        end
        attached2=1;
        xc(prot2)=ecmlocs(theNode2,1); % move the cortex points onto the fiber
        yc(prot2)=ecmlocs(theNode2,2);
    end
    if ((attached1 && attached2) || retracting) % adaptive timestep restriction
        dt=dt0/dtfactor;             % shorten up the timestep because the cell is stiffer
        if (t-tstar > largerdttime)  % after sufficient time put the timestep back to what it was
            dt=dt0;
        end
        if (~retracting)
            t0=0;
            maxprotlen=2*rc;
        end
    end
    % Force calculations. The nucleus / cortex force calculations are
    % different if the cell is attached
    if (attached1 && attached2)
        % Remesh the cortex and nucleus
        if (mod(t-tstar,0.01) < dt && t-tstar > dt)
            [newpos, x0,prot1,prot2,~,~] = remeshCortex(Nc,[xc yc],x0,prot1,prot2);
            [newnuc,nucref]=remeshNucleus(Nm,[xm ym],nucref);
            xc=newpos(:,1);
            yc=newpos(:,2);
            xm=newnuc(:,1);
            ym=newnuc(:,2);
        end
        % Elastic and bending forces are calculated differently
        fcel = calcnonunifelf(rc,xc,yc,x0(:,1),x0(:,2),gamc,gamc,kcsoft,kcrigid,prot1,prot2);
        fcel = fcel+calcnonunifbend(xc,yc,rc,x0,kbcrigid,prot1,prot2);
        % No bending force on the nucleus yet
        fmbend = zeros(Nm,2);
    else
        fcel = calcelasticforce(rc,xc,yc,gamc,kc);
        % If there is no nucleus bending rigidity, add a small bending
        % force when the cortex is making protrusions. 
        % Prevents nucleus from getting too deformed
        if (kbm < 1e-10)
            fmbend = calcbendingforce(xm,ym,rm,0.02);
        else
            fmbend = zeros(Nm,2);
        end
    end
    % Calculate elastic/bending forces with reference locations
    fmel=calcelasticforcerefs(rm,xm,ym,gamm,km,nucref); % nucleus elastic
    % Add the physical bending forces on the nucleus
    fmbend = fmbend+calcbendingforcerefs(xm,ym,rm,kbm,nucref);
    fnuctot=(fmel+fmbend)*hm; %Force/length
    fcortot=(fcel+pullingforce)*hc;
    % Mobility calculation
    alllocs=[xm ym; xc yc; ecmlocs];
    M=formMatrix(alllocs,eps,mu);
    if (rigidecm) % If the ECM is rigid, separate set of steps. Enforce the 0 velocity boundary 
                  % condition on the fibers, then solve the system for the
                  % motion of the cell.
        f=reshape([fnuctot; fcortot]',2*(Nc+Nm),1);
        wouldbeus=M(2*(Nc+Nm)+1:length(M),1:2*(Nc+Nm))*f;
        % Solve for the fiber force that keeps them rigid
        fibforce = M(2*(Nc+Nm)+1:length(M),2*(Nc+Nm)+1:length(M))\(-wouldbeus);
    else % otherwise the fiber force comes from elasticity
        fibforce=calcfibforce(ecmlocs,kteth,connections,ecmrefs);
    end
    totforce=[fnuctot; fcortot; fibforce];
    f=reshape(totforce',length(alllocs)*2,1);
    % Compute velocity
    uvals=reshape(M*f,size(alllocs'))';
    uvals=uvals+(log(bigR)-0.5)/(4*pi*mu)*sum(totforce); % add a constant velocity for the BC
    if (mod(floor(t/dt+1e-5),saveEvery)==0)
        plotmigcell; % Make a plot of what's happening
        % Calculate the area and aspect ratios of the cortex and nucleus
        AcorError = [AcorError calcArea(xc,yc) - pi*rc^2];
        AnucError = [AnucError calcArea(xm,ym) - pi*rm^2];
        CorAR = [CorAR (max(xc)-min(xc))/(max(yc)-min(yc))];
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
        nucpos=[nucpos; mean(xm) mean(ym)];
        [prot1, prot2] = pickprots(Nc,ppoints(:,length(jPtslist(1,:))+1),t0,km);
        jPtslist=[jPtslist [prot1; prot2]];
         % add the pinch and start over
        pullingforce=pullingforce+addpinch(fmag,Nc,sc,prot1)+addpinch(fmag,Nc,sc,prot2);  
    end
    % Take off the attachments if the cell has settled down. 
    if (max(sqrt(uvals(:,1).^2+uvals(:,2).^2)) < eps && attached1 && attached2)
        % move fiber
        [newpos, ~,prot1,prot2,np1,np2] = remeshCortex(Nc,[xc yc],x0,prot1,prot2);
        xc=newpos(:,1);
        yc=newpos(:,2);
        % Move the ECM fiber 2*eps away (like the cell "spitting it out")
        ecmlocs(theNode1,:)=ecmlocs(theNode1,:)+2*eps*np1;
        ecmlocs = correctPoints(ecmlocs,xc,yc,np1,2*eps);
        ecmlocs(theNode2,:)=ecmlocs(theNode2,:)+2*eps*np2;
        ecmlocs = correctPoints(ecmlocs,xc,yc,np2,2*eps);
        attached1=0;
        attached2=0;
        retracting=1;
        x0=x0_orig; % reset the elastic reference locations of the cortex
    end
    % Check if any ECM points are inside the cortex or cortex points inside
    % the nucleus and make the corrections by moving a distance epsilon
    corlocs = correctPoints([xc yc],xm,ym,[0 0],eps);
    xc = corlocs(:,1);
    yc = corlocs(:,2);
    ecmlocs = correctPoints(ecmlocs,xc,yc,[0 0],eps);
end
% Number of nucleus nodes ahead of ecm 32 and 39
% xline = (ecmlocs(39,1)-ecmlocs(32,1))/(ecmlocs(39,2)-ecmlocs(32,2))*(ym-ecmlocs(32,2))+ecmlocs(32,1);
% nucpen=sum(xm > xline)/Nm
% xline = (ecmlocs(39,1)-ecmlocs(32,1))/(ecmlocs(39,2)-ecmlocs(32,2))*(yc-ecmlocs(32,2))+ecmlocs(32,1);
% corpen=sum(xc > xline)/Nc



