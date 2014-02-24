% example of Thiem model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% define steady state
m.settime(0);

% set grid
m.setgrid(logspace(-1,3,41),1,true);

% set parameters
m.par.constant = [false(1,m.grid.nr-1), true];
m.par.kr = 50;

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
s = m.stress.q(1)/2/pi/m.par.kr * log(m.grid.r./m.grid.r(end));

% distance-drawdown graph
figure
semilogx(m.grid.r,m.s,'k-',m.grid.r,s,'kx')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
