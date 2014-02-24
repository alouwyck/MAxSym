% example of Butler model considering well with finite-thickness skin
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),1,true);

% set parameters
iR = 10;
m.par.kr = [ones(1,iR) 10*ones(1,m.grid.nr-iR)];
m.par.ss = [1e-5*ones(1,iR) 1e-3*ones(1,m.grid.nr-iR)];

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
it = 51:10:91;
[r,t] = meshgrid(m.grid.r,m.time.t(it));
R = m.grid.rb(iR+1);
b = r <= R;
s = zeros(size(r));
s(b) = m.stress.q(1)/2/pi * ...
       ((0.5772 + log(m.par.ss(end)*R*R/4/m.par.kr(end)./t(b)))...
        /2/m.par.kr(end) + log(r(b)/R)/m.par.kr(1));
b = ~b;
s(b) = m.stress.q(1)/4/pi/m.par.kr(end) * ...
       (0.5772 + log(m.par.ss(end)*r(b).^2/4/m.par.kr(end)./t(b)));
s(s>0) = 0;

% distance-drawdown graph
figure
semilogx(m.grid.r,squeeze(m.s(1,:,it))','-')
hold on
semilogx(m.grid.r,s','x')
set(gca,'fontsize',12)
xlabel('distance (m)')
ylabel('drawdown (m)')
xlim([1e-1 1e5])
ylim([-51 0])
legend(strcat(num2str(m.time.t(it),'%.0f'),'d'))
