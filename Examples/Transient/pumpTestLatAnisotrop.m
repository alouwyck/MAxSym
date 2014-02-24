% example of Hantush-Jacob model with laterally anisotropic conductivity
% A.LOUWYCK (2011)

% anisotropy
alfa = pi/6;
[kmax,kmin] = deal(20,5);
a = sqrt(kmax/kmin);
ke = sqrt(kmax*kmin);

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,4,91)));

% set grid
m.setgrid(logspace(-1,7,81),[1;1],true);

% set parameters
m.par.constant = [true; false];
m.par.kr = ke;
m.par.cz = 1000;
m.par.ss = 1e-3;

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2,1) = 100;

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% interpolate drawdown corresponding to apparent distances
r = 1;
theta = (30:30:120)/180*pi;
ra = r * sqrt(cos(alfa-theta).^2/a+sin(alfa-theta).^2*a);
sa = squeeze(m.interp1r(2,ra,[]))';

% analytical solution
ktheta = 1./(cos(alfa-theta).^2/kmax+sin(alfa-theta).^2/kmin);
[ktheta,t] = meshgrid(ktheta,m.time.t);
u = r*r./t * m.par.ss./ktheta/4;
v = r./sqrt(ktheta*m.par.cz);
for n = 1:numel(u)
    s(n) = -m.stress.q(2,1)/4/pi/ke * quadgk(@(y)exp(-y-v(n).*v(n)/4./y)./y, u(n), inf);
end
s = reshape(s,size(u));

% time-drawdown graph
figure
semilogx(m.time.t,sa,'-')
hold on
semilogx(m.time.t,s,'x')
set(gca,'fontsize',12)
title('r = 1m')
xlabel('time (d)')
ylabel('drawdown (m)')
legend(strcat(num2str(theta'/pi*180,'%.0f'),'°'))

% contour plot
[x,y] = meshgrid(-100:.5:100);
r = sqrt(x.^2+y.^2);
theta = atan2(y,x);
ra = r .* sqrt(cos(alfa-theta).^2/a+sin(alfa-theta).^2*a);
sa = reshape(m.interp1r(2,ra(:),length(m.time.t)),size(ra));
figure
contourf(x,y,sa);
axis square
colorbar
set(gca,'fontsize',12)
title(['t = ' num2str(m.time.t(end),'%.0f') 'd'])
xlabel('X (m)')
ylabel('Y (m)')
caxis([-9 0])