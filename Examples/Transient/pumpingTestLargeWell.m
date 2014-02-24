% example of pumping test in large-diameter well
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
rb = logspace(log10(0.5),7,801);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = 10;
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-3*ones(1,m.grid.nr-1)];

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model with wellbore storage
m.run;
sw = m.s;

% run model without wellbore storage
m.par.ss(1) = m.par.ss(2);
m.run;

% analytical solution
rw = m.grid.rb(2);
ir = 101:100:301;
r = [rw m.grid.r(ir)];

F = @(y) (y.*besselj(0,y)-2*m.par.ss(end)*besselj(1,y)).^2 + ...
         (y.*bessely(0,y)-2*m.par.ss(end)*bessely(1,y)).^2;

G = @(r,y) besselj(0,y*r/rw).*(y.*bessely(0,y)-2*m.par.ss(end)*bessely(1,y)) - ...
           bessely(0,y*r/rw).*(y.*besselj(0,y)-2*m.par.ss(end)*besselj(1,y));

H = @(r,y,t) (1-exp(-m.par.kr*t.*y.*y/rw^2/m.par.ss(end))) .* G(r,y)./y.^2./F(y);

for i = 1:length(r)
    for k = 1:length(m.time.t)
        s(i,k) = -2*m.stress.q(1)*m.par.ss(end)/pi/pi/m.par.kr * ...
                 quadgk(@(y)H(r(i),y,m.time.t(k)),0,Inf,'MaxIntervalCount',1500,'AbsTol',1e-5);
    end
end

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(sw(1,[1 ir],:))','-')
hold on
semilogx(m.time.t,squeeze(m.s(1,[1 ir],:))',':')
semilogx(m.time.t,s','x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(r','%.2f'),'m'))