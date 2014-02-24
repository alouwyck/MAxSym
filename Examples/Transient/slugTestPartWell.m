% example of slug test in partially penetrating well
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,0.1*ones(10,1),true);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:5,1) = true;
m.par.kr = 1;
m.par.kz = 0.1;
m.par.cz = [1e-8,zeros(1,m.grid.nr-1)];
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1)/5, 1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-6,500,5);

% run model
m.run;

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(6:10,1,:))')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
