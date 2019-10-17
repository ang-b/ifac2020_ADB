clearvars

A = [-2 -1; 1 0];
B = [1; 0];
C = ones(1, 2);
Ts = 0.1;
Aij = [0.1 0.3; 0 0];

nsys = 2;
testSys = lss();

link = {2, 1};

for i = 1:nsys
    testSys = testSys.addSystem(A,B,C,[], ['System ', int2str(i)]);
end

for i = 1:nsys
    for j = link{i}
        testSys = testSys.addCoupling(i,j,Aij);
    end
end

lsSim = @(t,x,u,d) full(testSys.getA)*x + full(testSys.getB)*u + d;
ssRange = @(x,v) sum(v(1:x-1)) + (1:v(x));

Tsim = 20;
t = 0:Ts:Tsim;
Nsamples = length(t);
ni = testSys.ni;

x0 = repmat([1;0], [nsys, 1]);
K = [0.5317, 2.5317];
L = K.';
x0hat = repmat([1;0], [nsys, 1]);

UIOs = cell(nsys, 1);
for i = 1:nsys
    UIOs{i} = UIO(A,B,C,Aij,L,Ts);
    % initial conditions default to 0 anyway
end

xsim = zeros(Nsamples, size(x0,1));
xsim(1,:) = x0.';
ysim = zeros(Nsamples, size(testSys.getC,1));
ysim(1,:) = testSys.getC*x0;
xhatsim = zeros(size(xsim));
xhatsim(1,:) = x0hat.';

for k=1:Nsamples-1
    
    u = zeros(nsys, 1);
    
    for i = 1:nsys
        si = testSys.getSystem(i);
        coupling = si.couplingA;
        nj = ni(coupling);
        xjstack = zeros(sum(nj),1);
        for j = 1:length(coupling)
            xjstack(ssRange(j,nj)) = x0hat(ssRange(coupling(j),ni));
        end
        u(i) = -K*x0hat(ssRange(i,ni)) - [1, 0]*full(si.Aij)*xjstack;
    end
        
   odesol = ode45(@(tt,x) lsSim(tt,x,u,0), [0 Ts], x0);
   
   x0 = odesol.y(:, end);
   xsim(k+1, :) = x0.';
   
   for i= 1:nsys
       x0hat(ssRange(i,ni)) = UIOs{i}.estimate(u(i), C*x0(ssRange(i,ni))).xhat;
   end
   xhatsim(k+1, :) = x0hat.';
end
 
% plot(t, xsim)
% hold on
% plot(t, xhatsim, 'LineStyle', '--');
% hold off


%% Test 1: estimation error converges to 0
assert(all(sum(abs(xhatsim(end:end-10,:) - xsim(end:end-10,:)),1) < eps*10));