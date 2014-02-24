% example of push-pull test model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(repmat(diff(logspace(-5,0,51)),[1 2]),50*ones(1,2));

% set grid
m.setgrid(logspace(-1,7,81),1,true);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-3;

% set stresses
m.stress(1).q = [-100, zeros(1,m.grid.nr-1)];
m.stress(2).q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% analytical solution
ir = 1:10:31;
[r,t] = meshgrid(m.grid.r(ir),m.time.t);
u = r.^2./t * m.par.ss/m.par.kr/4;
s1 = -1/4/pi/m.par.kr * expint(u);
s1(1) = 0;
ndt = [0; cumsum(m.time.ndt(:))];
dq = diff([0; sum(cat(1,m.stress.q),2)]);
s = zeros(size(s1));
for n = 1:length(ndt)-1
    s(ndt(n)+1:end,:) = s(ndt(n)+1:end,:) + dq(n)*s1(1:end-ndt(n),:);
end
 
% time-drawdown graph
figure
plot(m.time.t,squeeze(m.s(1,ir,:))','-')
hold on
plot(m.time.t,s,'x')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))
