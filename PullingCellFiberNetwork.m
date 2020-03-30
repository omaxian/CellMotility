% Simulation of a cell being pulled through the ECM matrix
% This corresponds to mechanism 1 in the paper
addpath('movecellfunctions');
close all;
% Initialize the user-defined parameters (separate file)
CellParameters_Mech1;
% Perform a bunch of operations to initialize all variables
dt=dt0;
dthc=2*pi/Nc;
dthm=2*pi/Nm;
kc=kcsoft;
kbc=0;
sc=[0:dthc:2*pi-dthc]';
sm=[0:dthm:2*pi-dthm]';
xm=rm*cos(sm);
ym=rm*sin(sm);
xc=rc*cos(sc);
yc=rc*sin(sc);
hc=dthc*rc;
hm=dthm*rm;
% Make the ECM matrix
if (loadecm)
    load(ecmfile); % If there is a file to load
else
    [fiberlocs,connections]=makeECM(numfibers,rc,rfib,nodebound,eps,fibspace); % clean, new ECM
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
attached=0;
reject=0;
% Select the point for the protrusion
% If ppoints(1) is -1, this will be random. Otherwise it 
% will just select the point in ppoints(1)
jPt = pickpoint(Nc,ppoints(1));
jPtslist=jPt;
nucpos=[0 0];
bigR=1000;
% Add the protrusion to the force vector pullingforce. 
pullingforce=pullingforce+addpinch(fmag,Nc,sc,jPt);
numProts=1;
retracting=0;
AcorError = 0;
AnucError = 0;
CorAR = 1;
NucAR = 1;
while (length(jPtslist) < length(ppoints))
    minnoded=10;
    % attached means if the cell is attached or not. So if it's not
    % attached, check how close the end of the protrusion is to a fiber.
    % retracting is if the protrusion went too far and it's coming back in
    if (~attached && ~retracting)
        [minnoded,theNode] = min(sqrt(sum(([xc(jPt) yc(jPt)]-ecmlocs).*([xc(jPt) yc(jPt)]-ecmlocs),2)));
        tstar=t;
        % Check how long the protrusion is getting 
        protlen=sqrt((xc(jPt)-mean(xc)).^2+(yc(jPt)-mean(yc)).^2)-rc;
    end 
    % Bring the protrusion back in if it gets too long (no ECM contact).
    if (protlen > maxprotlen)
        retracting=1;
        kc=kcrigid;
        pullingforce=pullingforce*0;
    end
    % If the min node distance is less than 2*epsilon (1 epsilon for the
    % fiber and 1 for the cortex node), bind the cortex node.
    if (minnoded < 2*eps || attached)
        attached=1;
        xc(jPt)=ecmlocs(theNode,1); % move the cortex points onto the fiber
        yc(jPt)=ecmlocs(theNode,2);
        pullingforce=pullingforce*0; % 0 out the pulling force
        kc=kcrigid;                  % increase the stiffness of the cortex.
        kbc=kbcrigid;                % add the bending rigidity to the cortex. 
    end
    if (retracting || attached) % timestep adjustment
        dt=dt0/dtfactor;
        if (t-tstar > largerdttime)  % after sufficient time put the timestep back to what it was
            dt=dt0;
        end
    end
    % Routine force calculations
    % Nucleus force
    fnucel=calcelasticforce(rm,xm,ym,gamm,km);
    if (~attached) % we give the nucleus a little bending rigidity when protrusions are forming
        fnucbend = calcbendingforce(xm,ym,rm,kbm);
        fnucel = fnucel+fnucbend;
    end
    fnuctot=fnucel*hm; %Force/length in 2D
    % Cortex forces
    fcel=calcelasticforce(rc,xc,yc,gamc,kc);
    fcbend = calcbendingforce(xc,yc,rc,kbc);
    fcortot=(fcel+fcbend+pullingforce)*hc; % Force/length in 2D
    % Stack locations
    alllocs=[xm ym; xc yc; ecmlocs];
    M=formMatrix(alllocs,eps,mu); % this is the most expensive routine
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
        kbc=0;
        nucpos=[nucpos; mean(xm) mean(ym)];
        jPt = pickpoint(Nc,ppoints(numProts+1)); % pick a new point and append to list
        numProts=numProts+1;
        jPtslist=[jPtslist jPt];
        pullingforce=pullingforce+addpinch(fmag,Nc,sc,jPt);   % add the pinch and start over
    end
    % Take off the attachment if the cell has settled down. 
    if (max(sqrt(uvals(:,1).^2+uvals(:,2).^2)) < eps && attached)
        % move fiber
        np1=getprotNormals(jPt,xc,yc);
        % Move the ECM fiber 2*eps away (like the cell "spitting it out")
        ecmlocs(theNode,:)=ecmlocs(theNode,:)+2*eps*np1/norm(np1);
        ecmlocs = correctPoints(ecmlocs,xc,yc,np1,2*eps);
        attached=0;
        retracting=1;
    end
    % Check if any ECM points are inside the cortex or cortex points inside
    % the nucleus and make the corrections by moving a distance epsilon
    ecmlocs = correctPoints(ecmlocs,xc,yc,[0 0],eps);
    corlocs = correctPoints([xc yc],xm,ym,[0 0],eps);
    xc = corlocs(:,1);
    yc = corlocs(:,2);
end