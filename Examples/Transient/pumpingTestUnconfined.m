% example of pumping test in unconfined layer
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),10,false);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-4;
m.par.sy = 0.2;

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,10,5);

% run model
m.run;

% analytical Theis solution
ir = 1:10:31;
[r,t] = meshgrid(m.grid.r(ir),m.time.t);
u = r.^2./t * m.par.sy/m.par.kr/m.grid.D/4;
s = -m.stress.q(1)/4/pi/m.par.kr/m.grid.D * expint(u);

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(1,ir,:))','-')
hold on
semilogx(m.time.t,s,':')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))
