% example of slug test with partially penetrating well in phreatic aquifer
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set solver
m.setsolver(1e-6,50,5);


% CONFINED CASE

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2),rb];
m.setgrid(rb,0.25*ones(4,1),true);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:2,1) = true;
m.par.kr = 1;
m.par.kz = 0.1;
m.par.cz = [1e-8,zeros(1,m.grid.nr-1)];
m.par.ss = [rb(2)^2*pi/m.grid.vol(1)/2,1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];

% run model
m.run;
s1 = squeeze(m.s(end,1,:));
disp(max(abs(m.totbud(:,1))))


% UNCONFINED CASE - CONSTANT HEAD

% set grid
m.setgrid(rb,[0.25;0.25*ones(4,1)],true);

% set parameters
m.par.constant = [true; false(4,1)];
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:3,1) = true;
m.par.kr = 1;
m.par.kz = 0.1;
m.par.cz = [1e-8,zeros(1,m.grid.nr-1)];
m.par.ss = [rb(2)^2*pi/m.grid.vol(1)/2,1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];

% run model
m.run;
s2 = squeeze(m.s(end,1,:));
disp(max(abs(m.totbud(:,1))))


% UNCONFINED CASE - SPECIFIC YIELD

% set grid
m.setgrid(rb,0.25*ones(4,1),false);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:2,1) = true;
m.par.kr = 1;
m.par.kz = 0.1;
m.par.cz = [1e-8,zeros(1,m.grid.nr-1)];
m.par.ss = [rb(2)^2*pi/m.grid.vol(1)/2,1e-5*ones(1,m.grid.nr-1)];
m.par.sy = 0.2;

% set stresses
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];

% run model
m.run;
s3 = squeeze(m.s(end,1,:));
disp(max(abs(m.totbud(:,1))))


% time-drawdown graph
figure
semilogx(m.time.t,[s1,s2,s3])
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
legend('confined','constant head','specific yield')
