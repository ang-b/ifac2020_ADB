%tstConvergence
clearvars, close all


%% Simulation setup
Tsim = 30;
Tsamp = 0.001;
t = 0:Tsamp:Tsim;
nTot= length(t);
Tref = 8;

xref = [0.1*sin(2*pi/Tref*t); zeros(1, nTot)];
% xref = zeros(2, nTot);

%% System setup -- continuous time: discretize subsys by subsys!
% Gravitational accelleration [m/s^2]
g = 9.80665;
% Length of inverted pendulum
l = .5;
% Number of subsystem matrices
N = 6;

% Edges in our LSS network
e = [1, 2;
     1, 3;
     2, 4;
     2, 6;
     3, 4;
     3, 5;
     4, 6;
     5, 6];

 % Number of edges
E = size(e,1);
% Masses
% m = [ .8 ;
%       .2 ;
%       .4 ;
%       .5 ];
m  = .1*ones(N,1);

%ceil(100*rand(N,1))/100;
% Height of interconnection
a = l*0.6;
% Interconnection 
kappa = [  27;
           40;
           35;
           53;
           33;
           42;
           25;
           38]*5;
% kappa = ceil(10*rand(E,1)) + 10;

% build LSS and static feedback control
Kpoles = [0.6 + 0.12i, 0.6 - 0.12i];

for i = N:-1:1
    subss(i) = buildSub(i,g,l,e,m,a,kappa,Tsamp);
    Kcell{i} = place(subss(i).A, subss(i).B, Kpoles);
%     Kcell{i} = [0.5/Tsamp*l/g 1];
end

K= blkdiag(Kcell{:});

LSSgraph = graph(e(:,1), e(:,2), kappa);

nI = size(subss(1).A, 1);
mI = size(subss(1).B, 2);
gI = size(subss(1).G, 2);
pI = size(subss(1).C, 1);

I = eye(nI*N);

% Full state dynamics
A = zeros(nI*N);
for i = 1:N
    for j = 1:N
        jj = find(subss(i).Ni == j);
        if size(jj,1) ~= 0
            A((1:nI)+nI*(i-1), (1:nI)+nI*(j-1)) = subss(i).Aij(:, (1:nI)+nI*(jj-1));
        elseif i == j 
            A((1:nI)+nI*(i-1), (1:nI)+nI*(i-1)) = subss(i).A;
        end
    end
end

B  = blkdiag(subss.B);
C  = blkdiag(subss.C);
G  = blkdiag(subss.G);
Ebar = blkdiag(subss.Ebar);

% edges between (system, controller)
Kedges = [ 1, 1 ; ...
           2, 2 ; ...
           3, 3 ; ...
           3, 4 ; ...
           4, 4 ; ...
           5, 5 ; ...
           6, 6];

%% Initialize state
% x = zeros(nI*N,nTot);
% y = zeros(pI*N,nTot);
x = zeros(nI, N, nTot);
y = zeros(pI, N, nTot);
u = zeros(mI, N, nTot);
utilde = zeros(mI, N, nTot);
yT = y;

% Initialized attack sequence
xA = zeros(nI, nTot);
yA = xA;

x0 = .3*rand(nI*N,1) - 0.15;
% set the initial speeds to zero (2nd component of state)
x0(repmat([false, true], 1, N)) = 0; 
u0 = K*x0;

x(:,:,1) = reshape(x0, nI, N, 1);
u(:,:,1) = reshape(u0, mI, N, 1);

%% attack initialization
mu = zeros(mI, nTot);
mui = mu;
muj = mu;
rho = mu;
nA = 3;
tA = 10;
kA = find(t == tA);
tA2 = 20;
kA2 = find(t == tA2);

nplantA =  length(Kedges(Kedges(:,1) == nA, 1));
% the system is unstable in open loop, therefore, to avoid numerical
% problems, the attacker designs her own controller

Ktildepoles = [0.7 + 0.2i, 0.7 - 0.2i];
Ktilde = place(subss(nA).A, subss(nA).B, Ktildepoles);

% this one destabilizes the system
% rho(:,kA:end) = 2*ones(1,mI*size(t(kA:end),1)).*(1-exp(-2*t(1:end-kA+1)));

% this one does not but induces bigger stationary oscillations
rho(:,kA:kA2) = 1*ones(mI,1).*sin(4*pi/Tref*t(kA:kA2)).*(1-exp(-2*t(1:kA2-kA+1)));


%% Observer initialization 

xd = zeros(nI,N,nTot);
xc = xd;
% xiEst = zeros(gI,N,nTot);

dcres = zeros(N,nTot);

yEst = xd;
%xd(:,:,1) = zeros(nI,N);
%yEst(:,:,1) = reshape(C*x0,pI,N)

% the observer dynamics is continuous, so place poles for CT
obsvPoles = [-50 -40];

for i = N:-1:1 % for local units
    % how many controllers in local unit i?
    ctrlrs_LU = Kedges(Kedges(:,2) == i, 1);
    nctrli = length(ctrlrs_LU);
    
    % subs p is controlled by LU(i).K{p}
    LU(i).K = cell(N, 1);
    for c=1:nctrli
        % how many plants are controller by controller c?
        nplantsc = length(Kedges(Kedges(:,1) == ctrlrs_LU(c), 1));
        LU(i).K{ctrlrs_LU(c)} = Kcell{ctrlrs_LU(c)} / nplantsc;
    end
    
    LU(i).UIO = UIO(subss(i).A, subss(i).B, subss(i).C, subss(i).Aij);
    LU(i).UIO.assignFPoles(obsvPoles);
    LU(i).UIO.tSpan = Tsamp;
    
%     Li = place(subss(i).A.', subss(i).C.', obsvPoles).';
    Li = LU(i).UIO.K1 + LU(i).UIO.H*subss(i).A;
    LU(i).LUE = CoupledLuenberger(subss(i).A, subss(i).B, subss(i).C, subss(i).Aij, Li, Tsamp);
end

%% Simulation

% cell i contains all controls to subss i
uj = cell(N, nTot);
utildej = uj;

for k = 1:nTot-1
    if mod(t(k+1), 1) == 0
        fprintf('Simulation time: % 6.2f/%.2f s\n', t(k+1), Tsim); 
    end
    % compute inputs
    
    mu(:,k) = -Ktilde * (xA(:,k) - rho(:,k));
    mui(:,k) = mu(:,k) / nplantA;
    muj(:,k) = mu(:,k) - mui(:,k);
    
    for i = 1:N % for subsystems  
        % how many controllers control plant i?
        ctrlrs_i = Kedges(Kedges(:,1) == i, 2);
        nctrli = length(ctrlrs_i);
        
        uol = zeros(mI, nctrli - 1);
        ulocal = zeros(mI, 1);
        
        uoltilde = uol;
        for c = 1:nctrli 
            if i == ctrlrs_i(c)
                u(:,i,k) = -LU(i).K{i} * (xd(:,i,k) - xref(:,k));
            else
                % ATTACK ON RECEIVED ESTIMATES GOES HERE
                % ctrlrs_i(c) receives xd(:,i,k) via communication network,
                % therefore an attack on the exchanged data alters such
                % estimate
                
                uol(:,c-1) = -LU(ctrlrs_i(c)).K{i} * (xd(:,i,k) - xref(:,k));
                uoltilde(:,c-1) = uol(:,c-1);
            end
                        
            % ATTACK ON RECEIVED INPUTS GOES HERE
            if i == nA && ctrlrs_i(c) ~= i
                uoltilde(:,c-1) = uol(:,c-1) + muj(:,k) ./ (nplantA - 1);
            end
        end       
        
        uj{i,k} = uol;
        utildej{i,k} = uoltilde;
    end
    
    utilde(:,:,k) = u(:,:,k);    
    utilde(:,nA,k) = u(:,nA,k) + mui(:,k);
    
    ujk = sum(uj{i,k}, 2);
    utildejk = sum(utildej{i,k}, 2);
    
    utildejk_stack = zeros(mI*N, 1);
    for i=1:N
         if isempty(utildej{i,k})
             utildejk_stack((1:mI)+mI*(i-1)) = zeros(mI, 1);
         else
             utildejk_stack((1:mI)+mI*(i-1)) = utildej{i,k};
         end
    end
    
    xA(:,k+1) = subss(nA).A*xA(:,k) + subss(nA).B*mu(:,k);
    yA(:,k) = subss(nA).C*xA(:,k);
    
    % System dynamics
    x(:,:,k+1) = reshape(...
                   A * reshape(x(:,:,k), [nI*N 1]) + ...
                   B * (reshape(utilde(:,:,k), [mI*N 1]) + utildejk_stack),  ... %finish reshape
                       [nI N 1]);
                   
    y(:,:,k) = reshape(C*reshape(x(:,:,k), [nI*N 1]), [pI N 1]);
    yT(:,:,k) = y(:,:,k);
    
    yT(:,nA,k) = y(:,nA,k) - yA(:,k);
  
    for i=1:N % for local units
        LU(i).UIO.estimate(u(:,i,k) + utildejk, yT(:,i,k));
        
        xd(:,i,k+1) = LU(i).UIO.xhat;
        
        xdj = reshape(xd(:,subss(i).Ni,k), [nI*length(subss(i).Ni) 1]);
        LU(i).LUE.estimate(u(:,i,k) + utildejk, yT(:,i,k), xdj);
        
        xc(:,i,k+1) = LU(i).LUE.xhat;
        
        dcres(i,k+1) = norm(xd(:,i,k+1) - xc(:,i,k+1));
    end
    
end