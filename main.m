%tstConvergence
clearvars, close all


%% Simulation setup
Tsim = 15;
Tsamp = 0.01;
t = 0:Tsamp:Tsim;
nTot= length(t);

xref = [0.1*sin(pi*t); zeros(1, nTot)];
% xref = zeros(2, nTot);

%% System setup -- continuous time: discretize subsys by subsys!
% Gravitational accelleration [m/s^2]
g = 9.80665;
% Length of inverted pendulum
l = .1;
% Number of subsystem matrices
N = 6;

% Edges in our LSS network
e = [1, 2;
     1, 4;
     2, 3;
     2, 5;
     3, 4;
     4, 5;
     4, 6];

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
           25]/10;
% kappa = ceil(10*rand(E,1)) + 10;

% build LSS and static feedback control
Kpoles = [0.7 + 0.12i, 0.7 - 0.12i];

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

Ai = blkdiag(subss.A);
B  = blkdiag(subss.B);
C  = blkdiag(subss.C);
G  = blkdiag(subss.G);
Ebar = blkdiag(subss.Ebar);

% edges between (system, controller)
Kedges = [ 1, 1 ; ...
           2, 2 ; ...
           3, 3 ; ...
           4, 4 ; ...
           4, 3 ; ...
           5, 5 ; ...
           6, 6];

%% Initialize state
% x = zeros(nI*N,nTot);
% y = zeros(pI*N,nTot);
x = zeros(nI, N, nTot);
y = zeros(pI, N, nTot);
u = zeros(mI, N, nTot);
yT = y;

% Initialized attack sequence
xA = zeros(nI, nTot);
xH = xA;
yA = xA;

x0 = 3*rand(nI*N,1) - 0.15;
% set the initial speeds to zero (2nd component of state)
x0(repmat([false, true], 1, N)) = 0; 
u0 = K*x0;

x(:,:,1) = reshape(x0, nI, N, 1);
u(:,:,1) = reshape(u0, mI, N, 1);

%% attack initialization
mu = zeros(gI, nTot);
nA = 3;
tA = 5;
kA = find(t == tA);

% this one destabilizes the system
mu(:,kA:end) = 0.05*ones(1,gI*size(t(kA:end),1)).*(1-exp(-2*t(1:end-kA+1)));

% this one does not but induces bigger stationary oscillations
% mu(:,kA:end) = .1*ones(gI,1).*sin(pi*t(kA:end)).*(1-exp(-2*t(1:end-kA+1)));

%% Observer initialization 

xd = zeros(nI,N,nTot);
xc = xd;
% xiEst = zeros(gI,N,nTot);

dcres = zeros(N,nTot);

yEst = xd;
%xd(:,:,1) = zeros(nI,N);
%yEst(:,:,1) = reshape(C*x0,pI,N)

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
    
    % design bank of UIOs
    LU(i).UIO = UIO(subss(i).A, subss(i).B, subss(i).C, subss(i).Aij);
    LU(i).UIO.assignFPoles(obsvPoles);
    LU(i).UIO.tSpan = Tsamp;
    
    Li = place(subss(i).A.', subss(i).C.', obsvPoles).';
    LU(i).LUE = CoupledLuenberger(subss(i).A, subss(i).B, subss(i).C, subss(i).Aij, Li, Tsamp);
end

%% Simulation

% cell i contains all controls to subss i
ucell = cell(N, nTot);

for k = 1:nTot-1
    if mod(t(k+1), 1) == 0
        fprintf('Simulation time: % 6.2f/%.2f s\n', t(k+1), Tsim); 
    end
    % compute inputs
    for i = 1:N % for subsystems
        % how many controllers control plant i?
        ctrlrs_i = Kedges(Kedges(:,1) == i, 2);
        nctrli = length(ctrlrs_i);
        
        uol = zeros(mI, nctrli);
        for c = 1:nctrli 
            % ATTACK ON RECEIVED ESTIMATES GOES HERE
            uol(:,c) = -LU(ctrlrs_i(c)).K{i} * (xd(:, i, k) - xref(:,k));  
        end
                
        % ATTACK ON RECEIVED INPUTS GOES HERE
        ucell{i,k} = uol;
        u(:,i,k) = sum(uol(:,c));
    end
    
    % Beginning of attack
    if k >= kA
        %Attack dynamics
        xA(:,k+1) = subss(nA).A*xA(:,k) + subss(nA).B*mu(:,k);
        yA(:,k) = subss(nA).C*xA(:,k) ;
        
        % System dynamics
        x(:,:,k+1) = reshape(...
                       A * reshape(x(:,:,k), [nI*N 1]) + ...
                       B * reshape(u(:,:,k), [mI*N 1]),  ... %finish reshape
                           [nI N 1]);
        x(:,nA,k+1) = x(:,nA,k) + subss(nA).B*mu(:,k);

        y(:,:,k) = reshape(C*reshape(x(:,:,k), [nI*N 1]), [pI N 1]);
        yT(:,:,k) = y(:,:,k);
        yT(:,nA,k) = y(:,nA,k) - yA(:,k);

    else
        % System dynamics
        x(:,:,k+1) = reshape( ...
                       A * reshape(x(:,:,k), [nI*N 1]) + ...
                       B * reshape(u(:,:,k), [mI*N 1]), ...
                           [nI N 1]);
        y(:,:,k) = reshape(C*reshape(x(:,:,k), [nI*N 1]), [pI N 1]);
        yT(:,:,k) = y(:,:,k); 
    end
    
    for i=1:N % for local units
        LU(i).UIO.estimate(u(:,i,k), yT(:,i,k));
        
        xd(:,i,k+1) = LU(i).UIO.xhat;
        
        xdj = reshape(xd(:,subss(i).Ni,k), [nI*length(subss(i).Ni) 1]);
        LU(i).LUE.estimate(u(:,i,k), yT(:,i,k), xdj);
        
        xc(:,i,k+1) = LU(i).LUE.xhat;
        
        dcres(i,k+1) = norm(xd(:,i,k+1) - xc(:,i,k+1));
    end
    
end