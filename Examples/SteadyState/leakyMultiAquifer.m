% example of steady-state Hemker model
% A.LOUWYCK (2011)

% input
T = [100; 10; 200];
c = [100; 500; 300; 1000];
Q = [100; 0; 1000];

% create Model object
m = MAxSym.Model;

% define steady state
m.settime(0);

% set grid
m.setgrid(logspace(-1,7,81),ones(5,1),true);

% set parameters
m.par.constant = [true; false(3,1); true];
m.par.kr = [0; T; 0];
m.par.cz = c;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2:end-1,1) = Q;

% set solver
m.setsolver(1e-5,100,5);

% run model
m.run;

% analytical solution
a = 1./T./c(1:end-1);
b = 1./T./c(2:end);
q = -Q/2/pi./T;
A = diag(a+b) - diag(a(2:end),-1) - diag(b(1:end-1),1);
[V,W] = eig(A);
for i = 1:m.grid.nr
    K = diag(besselk(0,m.grid.r(i)*sqrt(diag(W))));
    s(:,i) = V*K/V*q;
end

% distance-drawdown graph
figure
semilogx(m.grid.r,m.s(2:end-1,:)','-')
hold on
semilogx(m.grid.r,s,'x')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
legend(strcat('layer',num2str((1:3)')));