% example of combined injection and slug test
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
n = 300;
m.settime(repmat(diff(logspace(-5,0,n+1)),1,2),[n n]);

% set grid
rb = logspace(log10(0.03),6,801);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = 1;
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress(1).q = [-0.5, zeros(1,m.grid.nr-1)];
m.stress(2).q = m.stress(1).q;
m.stress(2).s0 = [1, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution (Theis)
rw = m.grid.rb(2);
t = m.time.t;
u = rw^2./t * m.par.ss(end)/m.par.kr/4;
s = -m.stress(1).q(1)/4/pi/m.par.kr * expint(u);

% analytical solution (Cooper et al.)
F = @(y) (y.*besselj(0,y) - 2*m.par.ss(end)*besselj(1,y)).^2 + ...
         (y.*bessely(0,y) - 2*m.par.ss(end)*bessely(1,y)).^2;
G = @(y,t) exp(-m.par.kr*t/rw^2.*y.*y/m.par.ss(end)) ./y./F(y);     
t = m.time.t(n+1:end)-m.time.t(n+1);
for i = 1:length(t)
    s(n+i) = s(n+i) + 8*m.stress(2).s0(1)*m.par.ss(end)/pi/pi * quadgk(@(y)G(y,t(i)),0,Inf);
end

% time-drawdown graph
figure
semilogx(t,squeeze(m.s(1,1,n+1:end))','k')
hold on
semilogx(t,s(n+1:end),'kx')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])