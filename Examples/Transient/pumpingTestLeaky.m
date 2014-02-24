% example of Hantush-Jacob model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),[1;1],true);

% set parameters
m.par.constant = [true; false];
m.par.kr = 10;
m.par.cz = 1000;
m.par.ss = 1e-3;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2,1) = 100;

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
ir = 1:10:31;
[r,t] = meshgrid(m.grid.r(ir),m.time.t);
u = r.*r./t * m.par.ss/m.par.kr/4;
for n = 1:numel(u)
    v = r(n)/sqrt(m.par.kr*m.par.cz);
    s(n) = -m.stress.q(2,1)/4/pi/m.par.kr * quadgk(@(y)exp(-y-v.*v/4./y)./y, u(n), inf);
end
s = reshape(s,size(u));

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(2,ir,:))','-')
hold on
semilogx(m.time.t,s,'x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))

% plot of leakage and storage fraction vs time
figure
semilogx(m.time.t(2:end),[sum(squeeze(m.qs(2,:,:))); sum(squeeze(m.qz(2,:,2:end)))]'/-m.stress.q(2,1))
set(gca,'fontsize',12)
xlabel('time (d)')
ylim([0 1])

