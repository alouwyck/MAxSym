% example of step-drawdown test considering non-linear well loss
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(repmat(diff(logspace(-5,0,51)),[1 4]),50*ones(1,4));

% set grid
rb = logspace(-1,7,801);
m.setgrid([rb(1)^2/rb(2), rb],1,true);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-3;

% set stresses
m.stress(1).q = sparse(1,1,100,1,m.grid.nr);
m.stress(2).q = sparse(1,1,200,1,m.grid.nr);
m.stress(3).q = sparse(1,1,300,1,m.grid.nr);
m.stress(4).q = sparse(1,1,400,1,m.grid.nr);

% set solver
m.setsolver(1e-5,10);

% run model
m.run;

% add well loss
sw = squeeze(m.s(1,1,:));
ndt = [0; cumsum(m.time.ndt(:))];
dq = diff([0; sum(cat(1,m.stress.q),2)]);
for n = 1:length(ndt)-1
    sw(ndt(n)+2:ndt(n+1)+1) = sw(ndt(n)+2:ndt(n+1)+1) - 1e-5*m.stress(n).q(1)^2;
end

% time-drawdown graph
figure
plot(m.time.t,squeeze(m.s(1,1,:))','k-',m.time.t,sw,'r-')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend('not corrected','corrected')
