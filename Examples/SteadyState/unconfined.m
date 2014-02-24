% example of Thiem-Dupuit model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% define steady state
m.settime(0);

% set grid
m.setgrid(logspace(-1,3,41),10,false);

% set parameters
m.par.constant = [false(1,m.grid.nr-1), true];
m.par.kr = 5;

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,100,5);

% run model
m.run;

% analytical solution
s = sqrt(m.grid.D^2 - m.stress.q(1)/pi/m.par.kr * log(m.grid.r(end)./m.grid.r)) - m.grid.D;

% distance-drawdown graph
figure
semilogx(m.grid.r,m.s,'k-',m.grid.r,s,'kx')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
