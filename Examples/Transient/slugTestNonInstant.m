% example of slug test with non-instantaneous head change
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
n = 31;
ti = linspace(0,3,n)/24/60/60; % 3 seconds
m.settime([diff(ti),diff(logspace(-6,1,251))],[ones(1,n-1),250]);

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = 1;
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-5*ones(1,m.grid.nr-1)];

% set stresses
for k = 1:n-1
    m.stress(k).s0 = [1/(n-1) ,zeros(1,m.grid.nr-1)];
end

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
F = @(y) (y.*besselj(0,y) - 2*m.par.ss(end)*besselj(1,y)).^2 + ...
         (y.*bessely(0,y) - 2*m.par.ss(end)*bessely(1,y)).^2;
G = @(y,t) exp(-m.par.kr*t/m.grid.rb(2)^2.*y.*y/m.par.ss(end)) ./y./F(y);
t0 = m.time.t(n);
for k = n:length(m.time.t)
    s(k-n+1) = 8*m.par.ss(end)/pi/pi * quadgk(@(y)G(y,m.time.t(k)-t0),0,Inf);
end

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(1,1,:)),'k',m.time.t(n:end),s,'r')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
