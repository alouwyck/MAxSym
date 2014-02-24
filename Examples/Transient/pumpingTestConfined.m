% example of Theis model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),1,true);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-3;

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
ir = 1:10:31;
[r,t] = meshgrid(m.grid.r(ir),m.time.t);
u = r.^2./t * m.par.ss/m.par.kr/4;
s = -m.stress.q(1)/4/pi/m.par.kr * expint(u);

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(1,ir,:))','-')
hold on
semilogx(m.time.t,s,'x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))

% interpolate drawdown
r = [1 5 10];
t = 1:100;
si = m.interp2(1,r,t);

% analytical solution
[t,r] = meshgrid(t,r);
u = r.^2./t * m.par.ss/m.par.kr/4;
s = -m.stress.q(1)/4/pi/m.par.kr * expint(u);

% plot interpolated drawdown vs time
figure
plot(t',squeeze(si)','k-',t',s','kx')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')

