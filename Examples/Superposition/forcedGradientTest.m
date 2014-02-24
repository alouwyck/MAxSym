% example of forced gradient test model
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,1,61)));

% set grid
m.setgrid(logspace(-1,7,81),10,true);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-4;

% set stresses
m.stress.q = [1, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% time-drawdown graph for extraction well
figure
sw = 100 * (m.s(1,1,:) - m.interp1r(1,50,[]));
semilogx(m.time.t,squeeze(sw))
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')

% 2D raster and well positions
[x,y] = meshgrid(1:100);
[xp,yp] = deal(50,25);
[xi,yi] = deal(50,75);

% interpolation
rp = sqrt((x-xp).^2+(y-yp).^2);
sp = squeeze(m.interp1r(1,rp(:),length(m.time.t)));
ri = sqrt((x-xi).^2+(y-yi).^2);
si = squeeze(m.interp1r(1,ri(:),length(m.time.t)));

% superposition
s = reshape(100*(sp-si),size(rp));

% contour plot
figure
contourf(x,y,s)
colorbar
axis equal
xlim([1 100])
ylim([1 100])
caxis([-1 1])
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('y (m)')

% effective velocities
n = 0.3;
[x,y] = deal(x(2:end-1,2:end-1),y(2:end-1,2:end-1));
qx = -m.par.kr/n*diff(s,1,2);
qy = -m.par.kr/n*diff(s,1,1);
qx = (qx(2:end-1,1:end-1) + qx(2:end-1,2:end))/2;
qy = (qy(1:end-1,2:end-1) + qy(2:end,2:end-1))/2;
ip = find(x==xp & y==yp);
ii = find(x==xi & y==yi);
di = [0 1 -1 98 -98];
qx([ip+di ii+di]) = 0;
qy([ip+di ii+di]) = 0;

% quiver plot
figure
for i = 1:2
    subplot(1,2,i)
    quiver(x,y,qx/2,qy/2,0)
    axis equal
    xlim([40 60])
    ylim([15 35]+(i-1)*50)
    set(gca,'fontsize',12)
    xlabel('x (m)')
    ylabel('y (m)')
end