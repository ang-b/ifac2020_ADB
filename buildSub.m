function [ subss ] = buildSub ( i , g , l , e , m , a , k , tSamp )
% Build subsystem i, according to continuous time equation of 
% interconnected inverted pendulums (1.61) in Siljak, Decentralized control
% of Complex Systems.
%
% Return subsystem i
%
% Parameters:
%       - i     : index of subsystem;
%       - g     : gravitational accelleration;
%       - l(i)  : length of pendulum i;
%       - e     : set of all edges (undirected);
%       - m(i)  : mass of subsystem i;
%       - a     : height of interconnecting spring;
%       - k(E)  : spring constant of interconnection index E;

A = [   0     1 ; 
        -g/l  0 ];

I = eye(2);

% A = [  0   1 ; 
%       -2  -3 ];

B = [ 0 ; -1/(m(i)*l^2) ];
% B = [ 0;1 ];
mI = size(B,2);

G =  [ 0;1 ];
gI = size(G,2);

C = eye(2);
    
Ei1 = find(e(:,1) == i); % find indices of outbound neighbors
Ei2 = find(e(:,2) == i); % find indices of inbound neighbors
Ei  = [Ei2;Ei1];

Ni  = [e(Ei2,1);e(Ei1,2)];
NNi = length(Ni);
nIJ = 2*NNi;

Aij = zeros(2,2*NNi);

Aij(2,([1:NNi]-1)*2+1) = k(e(Ei))*a^2 / (m(i)*l^2);

A = A - [sum(Aij,2),zeros(2,1)];

% Btmp = [ B , Aij ];

% tmp = ss(A,Btmp,C,[]);
% ssD = c2d(tmp,tSamp);

% Euler approximation
ssD.A   = (I + tSamp*A);
ssD.B   = tSamp*B;
ssD.Aij = tSamp*Aij;

subss.A    = ssD.A;
subss.B    = ssD.B;
subss.Aij  = ssD.Aij;
subss.Ebar = ssD.Aij(2,:);
% subss.B    = ssD.B(:,1:mI);
% subss.Aij  = ssD.B(:,mI+1:mI+nIJ);
% subss.Ebar = ssD.B(2,mI+1:mI+nIJ);
subss.G    = G;
subss.C    = C;
subss.Ni   = Ni;

if subss.G*subss.Ebar ~= subss.Aij
    error('You have done something wrong')
end

