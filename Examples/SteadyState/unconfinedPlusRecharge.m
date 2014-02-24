% example of Thiem-Dupuit model with recharge
% A.LOUWYCK (2011)

% radius of influence
Q = 100;
N = 5e-4;
R = sqrt(Q/pi/N);

% create Model object
m = MAxSym.Model;

% define steady state
m.settime(0);

% set grid
rb = logspace(-1,log10(R),31);
m.setgrid([rb,rb(end)+1e-5],10,false);

% set parameters
m.par.constant = [false(1,m.grid.nr-1), true];
m.par.kr = 5;

% set stresses
m.stress.q = [Q, -N*m.grid.hs(2:end-1), 0];

% set solver
m.setsolver(1e-5,100,5);

% run model
m.run;

% analytical solution
s = sqrt(m.grid.D^2 - Q/2/pi/m.par.kr * (log(Q/pi/N./(m.grid.r).^2)-1) - N/2/m.par.kr*m.grid.r.^2) - m.grid.D;

% distance-drawdown graph
figure
semilogx(m.grid.r,m.s,'k-',m.grid.r,s,'kx')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
