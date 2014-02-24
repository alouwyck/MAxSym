% example of two slug tests performed in a confined aquifer
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2),rb];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = 1;
m.par.ss = [rb(2)^2*pi/m.grid.vol(1),1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% drawdown within well
sw = squeeze(m.s(1,1,:));

% drawdown at 5 m
sr = squeeze(m.interp1r(1,5,[]));

% initial head changes
H0 = [1.5; 2];

% superposition
s = [sw,sr] * [H0,flipud(H0)];

% time-drawdown plot
figure
semilogx(m.time.t,s)
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
legend('well 1','well 2')
