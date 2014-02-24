% example of steady-state multi-aquifer model with recharge
% A.LOUWYCK (2011)

% input
T = [50; 100; 10];
c = [100; 250];
Q = [0; 500; 0];
N = 5e-4;
Rc = sqrt(sum(Q)/pi/N);
 
% create Model object
m = MAxSym.Model;
 
% define steady state
m.settime(0);
 
% set grid
rb1 = logspace(-1,log10(Rc),41);
rb2 = logspace(log10(Rc),7,41);
m.setgrid([rb1,rb2(2:end)],ones(3,1),true);
 
% set parameters
m.par.constant = false(m.grid.nz,m.grid.nr);
m.par.constant(1,end) = true;
m.par.kr = T;
m.par.cz = c;
 
% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(:,1) = Q;
m.stress.q(1,1:40) = -N*m.grid.hs(1:40);
 
% set solver
m.setsolver(1e-5,100,5);
 
% run model
m.run;
 
% analytical solution
c = [Inf; c(:); Inf];
a = 1./T./c(1:end-1);
b = 1./T./c(2:end);
A = diag(a+b) - diag(a(2:end),1) - diag(b(1:end-1),-1);
[V,W] = eig(A);
W = diag(W);
[~,i] = min(W);
W(i) = [];
W = repmat(sqrt(W(:)'),3,1);
V(:,i) = [];
t = T/sum(T);
 
b = true(3,1);
[Q,i] = max(Q);
b(i) = false;
X = V(b,:)\(Q*t(b));
X = repmat(X',3,1);
Y = (V(2:end,:)/Rc.*W(2:end,:))\(-N*t(2:end));
Y = repmat(Y',3,1);
 
for i = 1:m.grid.nr
    x = X/2/pi .* besselk(0,m.grid.r(i)*W) .* V;
    sw(:,i) = (Q/2/pi * log(m.grid.r(i)) * t + sum(x,2)) ./ T;
    if m.grid.r(i) <= Rc
        y = Y.*V .* (1./W/Rc - ... 
		     besselk(1,Rc*W).*besseli(0,m.grid.r(i)*W));
        sc(:,i) = (-N/4 * (m.grid.r(i)^2-Rc^2)*t + sum(y,2)) ./ T;
    else
        y = Y.*V .* besseli(1,Rc*W).*besselk(0,m.grid.r(i)*W);
        sc(:,i) = (-N*Rc^2/2*log(m.grid.r(i)/Rc)*t + sum(y,2)) ./ T;
    end
end
 
s = sw + sc -(sw(1,end)+sc(1,end));
 
% distance-drawdown graph
figure
semilogx(m.grid.r,m.s','-')
hold on
semilogx(m.grid.r,s','x')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
legend(strcat('layer',num2str((1:3)')));
