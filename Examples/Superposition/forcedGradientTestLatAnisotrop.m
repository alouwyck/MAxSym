% example of forced gradient test conducted in laterally anisotropic aquifer
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

% 2D raster
[x,y] = meshgrid(1:100);

% wells
[xp,yp] = deal(50,25);
[xi,yi] = deal(50,75);

% interpolation
alfa = pi/6;
a = sqrt(20/5);

r = sqrt((x-xp).^2+(y-yp).^2);
theta = atan2(y-yp,x-xp);
r = r .* sqrt(cos(alfa-theta).^2/a+sin(alfa-theta).^2*a);
sp = m.interp1r(1,r(:),length(m.time.t));
r = sqrt((x-xi).^2+(y-yi).^2);
theta = atan2(y-yi,x-xi);
r = r .* sqrt(cos(alfa-theta).^2/a+sin(alfa-theta).^2*a);
si = m.interp1r(1,r(:),length(m.time.t));

% superposition
s = reshape(100*(sp-si),size(r));

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

kx = 1/(cos(alfa)^2/20+sin(alfa)^2/5);
qx = -kx/n*diff(s,1,2);
qx = (qx(2:end-1,1:end-1) + qx(2:end-1,2:end))/2;

ky = 1/(cos(alfa-pi/2)^2/20+sin(alfa-pi/2)^2/5);
qy = -ky/n*diff(s,1,1);
qy = (qy(1:end-1,2:end-1) + qy(2:end,2:end-1))/2;

ip = find(x==xp & y==yp);
ii = find(x==xi & y==yi);
di = [0 1 -1 98 -98];
qx([ip+di ii+di]) = 0;
qy([ip+di ii+di]) = 0;

% quiver plot around wells
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

% quiver plot entire grid
[ip,jp] = ind2sub(size(x),ip);
[ii,ji] = ind2sub(size(x),ii);
di = -5:5;
qx(ip+di,jp+di) = 0;
qx(ii+di,ji+di) = 0;
qy(ip+di,jp+di) = 0;
qy(ii+di,ji+di) = 0;

figure
i = 1:2:size(x,1);
quiver(x(i,i),y(i,i),qx(i,i),qy(i,i),3)
axis equal
xlim([0 100])
ylim([0 100])
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('y (m)')
