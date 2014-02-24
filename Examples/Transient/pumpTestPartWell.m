% example of Hantush-Weeks model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),.1*ones(10,1),true);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:5,1) = true;
m.par.kr = 10;
m.par.kz = 1;
m.par.cz = sparse(m.grid.nz-1,m.grid.nr);
m.par.cz(6:end,1) = 1e-5;
m.par.ss = 1e-3;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(6:end,1) = 20;

% set solver
m.setsolver(1e-5,1000,5);

% run model
m.run;

% analytical solution
ir = 11:10:31;
ilay = [3 8];
[t,r] = meshgrid(m.time.t,m.grid.r(ir));
u = r.*r./t * m.par.ss/m.par.kr/4;
D = sum(m.grid.D);
[b,d] = deal(D,D/2);
[bq,dq] = deal(ilay/m.grid.nz*D,(ilay-1)/m.grid.nz*D);
a = sqrt(m.par.kz/m.par.kr);

for i = 1:length(ilay)
    F = zeros(size(u)); 
    for j = 1:numel(u)
        for n = 1:5
            v = n*pi*a*r(j)/D;
            F(j) = F(j) + quadgk(@(y)exp(-y-v.*v/4./y)./y, u(j),inf)/n/n * ...
                   (sin(n*pi*b/D)-sin(n*pi*d/D))*(sin(n*pi*bq(i)/D)-sin(n*pi*dq(i)/D));
        end
    end
    s(i,:,:) = -sum(m.stress.q(:,1))/4/pi/m.par.kr/D * ...
               (expint(u) + 2*D*D/pi/pi/(b-d)/(bq(i)-dq(i)) * F);
end

% time-drawdown graph
figure
for i = 1:length(ilay)
    subplot(1,2,i)
    semilogx(m.time.t,reshape(m.s(ilay(i),ir,:),length(ir),[])','-')
    hold on
    semilogx(m.time.t,reshape(s(i,:,:),length(ir),[])','x')
    set(gca,'fontsize',12)
    title(['layer ' int2str(ilay(i))])
    xlabel('time (d)')
    ylabel('drawdown (m)')
    ylim([-16 0])
    legend(strcat(num2str(m.grid.r(ir)','%.2f'),'m'))
end

% analytical solution for drawdown in the pumping well
r = m.grid.r(1);
u = r*r./m.time.t * m.par.ss/m.par.kr/4;
F = zeros(size(u)); 

for j = 1:numel(u)
    for n = 1:10
        v = n*pi*a*r/D;
        F(j) = F(j) + quadgk(@(y)exp(-y-v.*v/4./y)./y, u(j),inf)/n/n * ...
               (sin(n*pi*b/D)-sin(n*pi*d/D))^2;
    end
end

s = -sum(m.stress.q(:,1))/4/pi/m.par.kr/D * ...
   (expint(u) + 2*D*D/pi/pi/(b-d)^2 * F);

% time-drawdown graph
figure
semilogx(m.time.t,squeeze(m.s(6:end,1,:))','k',m.time.t,s,'kx')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
