% example of step-drawdown test in well with infinitesimal skin
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
m.par.cr = sparse(1,1,1e-2,1,m.grid.nr-1);

% set stresses
m.stress(1).q = sparse(1,1,100,1,m.grid.nr);
m.stress(2).q = sparse(1,1,200,1,m.grid.nr);
m.stress(3).q = sparse(1,1,300,1,m.grid.nr);
m.stress(4).q = sparse(1,1,400,1,m.grid.nr);

% set solver
m.setsolver(1e-5,10);

% run model
m.run;

% analytical solution
r = [m.grid.rb(2), m.grid.r(101:100:301)];
[r,t] = meshgrid(r,m.time.t);
u = r.^2./t * m.par.ss/m.par.kr/4;
s1 = -1/4/pi/m.par.kr * expint(u);
s1(1) = 0;
ndt = [0; cumsum(m.time.ndt(:))];
dq = diff([0; sum(cat(1,m.stress.q),2)]);
s = zeros(size(s1));
for n = 1:length(ndt)-1
    s(ndt(n)+1:end,:) = s(ndt(n)+1:end,:) + dq(n)*s1(1:end-ndt(n),:);
    ds(:,n) = m.par.cr(1)*m.stress(n).q(1)/2/pi/m.grid.rb(2);
    s(ndt(n)+2:ndt(n+1)+1,1) = s(ndt(n)+2:ndt(n+1)+1,1) - ds(:,n);
end

% time-drawdown graph
figure
plot(m.time.t,squeeze(m.s(1,1:100:301,:))','-')
hold on
plot(m.time.t,s,'x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(r(1,:)','%.2f'),'m'))

