% example of pumping test in unconfined layer with recharge
% A.LOUWYCK (2011)

% radius of influence
Q = 100;
N = 5e-4;
R = sqrt(Q/pi/N);

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
rb = logspace(-1,log10(R),31);
m.setgrid([rb(1:end-1),logspace(log10(R),7,51)],10,false);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-4;
m.par.sy = 0.2;

% set stresses
m.stress.q = [Q, -N*m.grid.hs(2:30), zeros(1,50)];

% set solver
m.setsolver(1e-5,10,5);

% run model
m.run;

% time-drawdown graph
ir = 1:10:21;
figure
semilogx(m.time.t,squeeze(m.s(1,ir,:))','-')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))

% distance-drawdown graph
it = 31:10:91;
figure
semilogx(m.grid.r,squeeze(m.s(1,:,it))','-')
hold on
semilogx([R R],ylim(gca),'k:')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.time.t(it),'%.0e'),'d'))
xlim([0.1 1e3])

% plot of storage fraction vs time
figure
qs = squeeze(m.qs);
b = qs >= 0;
qs(b) = 0;
semilogx(m.time.t(2:end),-sum(qs)/Q)
set(gca,'fontsize',12)
xlabel('time (d)')
xlim([1e-3 1e4])
ylim([0 1])