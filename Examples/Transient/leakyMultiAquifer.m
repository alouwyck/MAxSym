% example of transient Hemker model
% A.LOUWYCK (2011)

% input
T = [100; 10; 200];
c = [100; 500; 300; 1000];
Q = [100; 0; 1000];
d = 1e4;
S = T/d;

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),ones(5,1),true);

% set parameters
m.par.constant = [true; false(3,1); true];
m.par.kr = [0; T; 0];
m.par.ss = [0; S; 0];
m.par.cz = c;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2:end-1,1) = Q;

% set solver
m.setsolver(1e-8,10,5);

% run model
m.run;

% analytical solution
ir = [11 31];
a = 1./S./c(1:end-1);
b = 1./S./c(2:end);
q = -Q/2/pi./T;
A = diag(a+b) - diag(a(2:end),-1) - diag(b(1:end-1),1);
[X,Y] = eig(A);
Y = diag(Y);
for i = 1:length(m.time.t)
    t = m.time.t(i);
    for j = 1:length(ir)
        r = m.grid.r(ir(j));
        for n = 1:length(Y)
            v = r*sqrt(Y(n)/d);
            W(n) = quadgk(@(y)exp(-y-v.*v/4./y)./y,r*r/4/t/d,Inf)/2;
        end
        s(:,j,i) = X*diag(W)/X*q;
    end
end

% time-drawdown graphs
figure
for i = 1:length(ir)
    subplot(1,2,i)
    semilogx(m.time.t,squeeze(m.s(2:end-1,ir(i),:))','-')
    hold on
    semilogx(m.time.t,squeeze(s(:,i,:))','x')
    set(gca,'fontsize',12)
    title(['r = ' num2str(m.grid.r(ir(i)),'%.2f') 'm'])
    xlabel('time (d)')
    ylabel('drawdown (m)')
    legend(strcat('layer',num2str((1:3)')));
end

% plot of leakage and storage fraction vs time for each aquifer
f = squeeze(sum(m.bud([1 5],:,:),2))';
f(:,3:5) = squeeze(-sum(m.qs(2:4,:,:),2))';
f = f/sum(Q(:,1));
figure
semilogx(m.time.t(2:end),f);
set(gca,'fontsize',12)
xlabel('time (d)')
legend({'leakage top';'leakage bottom';'storage aquifer 1';...
        'storage aquifer 2';'storage aquifer 3';});
