% example of De Glee model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% define steady state
m.settime(0);

% set grid
m.setgrid(logspace(-1,7,81),[1;1],true);

% set parameters
m.par.constant = [true;false];
m.par.kr = [0;10];
m.par.cz = 1000;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2,1) = 100;

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
s = -m.stress.q(2,1)/2/pi/m.par.kr(2) * besselk(0,m.grid.r/sqrt(m.par.kr(2)*m.par.cz));

% distance-drawdown graph
figure
semilogx(m.grid.r,m.s(2,:),'k-',m.grid.r,s,'kx')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
