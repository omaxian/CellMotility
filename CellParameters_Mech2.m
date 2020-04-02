% The parameters (user-defined) for the moving cell through ECM
% simulations. These are the parameters for mechanism 2. 
% Fluid parameters
mu=1;               % fluid viscosity
eps=0.075;          % for regularized Stokeslets and also how close cell gets to fiber 
% Parameters for the cell (cortex & membrane)
Nc=120;             % the number of points for the cortex
Nm=80;              % number of nuclear points
gamm=0;             % surface tension term for the nucleus
km=Kmvals(iRegime); % nuclear stiffness
kbm = 0.0;          % physical bending stiffness of the nucleus
kbcrigid = 0.0;     % bending stiffness of cortex
gamc=0;             % surface tension term for the cortex
kcsoft=1;           % stiffness of the cortex when it is loose and floppy
kcrigid=Kcvals(iRegime);        % stiffness of the cortex when it stiffens up
rc=0.5;             % radius of the cortex
rm=0.45;            % radius of the nucleus
% Parameters related to the protrusions
fmag=500;           % magnitude of the protrusive force
maxprotlen=3*rc;    % initial maximum length of a protrusion (2*rc for mechanism 2)
%ppoints = [92 105 -1; 112 6 -1]; % for sparse
% ppoints = [115 114 -1; 13 10 -1]; % for Figs. S10-S16...
ppoints = -ones(2,protsperTrial+1); % for random trials
% Parameters for the ECM
fibspace=0.5*rc;    % minimum spacing of the fiber centers (need 2*rfib so they don't overlap)
numfibers=60;       % number of fibers in the ECM
rigidecm=0;         % whether the ECM is rigid or not
kteth=50;           % stiffness of the springs that make up the ECM lattice
nodebound=8*rc;     % the box where the ECM nodes can go
loadecm=1;          % if we are loading an existing ECM
ecmfile='ecm60.mat';% file name of the existing ECM
% Parameters related to timestepping
dt0=0.001;          % the base timestep
dtfactor=5;         % the factor the timestep goes down by when the cortex gets stiff
saveEvery=100;       % how often to save the frame for a movie
largerdttime=0.05;  % how long the shorter timestep is needed
