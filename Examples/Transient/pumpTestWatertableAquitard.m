% example of Bouton-Cooley model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),[2;1],false);

% set parameters
m.par.kr = [0.1;100];
m.par.cz = 100;
m.par.ss = [1e-5;1e-3];
m.par.sy = 0.2;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2,1) = 100;

% set solver
m.setsolver(1e-5,10,5);

% run model
m.run;

% Boulton solution
ir = 1:10:31;
alfa = 1/m.par.sy/m.par.cz;
L = sqrt(m.par.kr(2)*m.par.cz);
eta = 1 + m.par.sy/m.par.ss(2);
for i = 1:length(ir)
    r = m.grid.r(ir(i));
    for k = 1:length(m.time.t)
        t = m.time.t(k);
        F = @(x) besselj(0,r.*x/L)./x .* ...
                 (1-1./(x.^2+1).*exp(-alfa*t*x.^2./(x.^2+1))-...
                  x.^2./(x.^2+1).*exp(-alfa*t*eta*(x.^2+1)));
        s(k,i) = -m.stress.q(2,1)/2/pi/m.par.kr(2) * ...
                 quadgk(F,0,Inf,'AbsTol',5e-3);
    end
end

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(2,ir,:))','-')
hold on
semilogx(m.time.t,s,'x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))
