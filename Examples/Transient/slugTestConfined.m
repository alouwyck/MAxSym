% example of Cooper et al. model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = 1;
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1 ,zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
F = @(y) (y.*besselj(0,y) - 2*m.par.ss(end)*besselj(1,y)).^2 + ...
         (y.*bessely(0,y) - 2*m.par.ss(end)*bessely(1,y)).^2;
G = @(y,t) exp(-m.par.kr*t/m.grid.rb(2)^2.*y.*y/m.par.ss(end)) ./y./F(y);     
for n = 1:length(m.time.t)
    s(n) = 8*m.stress.s0(1)*m.par.ss(end)/pi/pi * quadgk(@(y)G(y,m.time.t(n)),0,Inf);
end

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(1,1,:)),'k',m.time.t,s,'kx')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
