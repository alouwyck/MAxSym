% example of model for multi-level extraction near lateral boundaries
% A.LOUWYCK (2011)


% MAP

% coordinates and sign for well and image wells
[xp,yp,q] = deal([85 85 115 115],[85 115 85 115],[1 1 -1 -1]);

% map indicating well, image wells and lateral bounds
figure
plot([70 100; 100 100]',[100 100; 70 100]');
axis equal
xlim([70 130])
ylim([70 130])
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('y (m)')

hold on
for i = 1:4
    if i==1
        col = 'b';
    else
        col = 'r';
    end
    if q(i)>0
        sgn = '';
    else
        sgn = '-';
    end
    plot(xp(i),yp(i),strcat(col,'o'))
    text(xp(i)+1,yp(i),strcat(sgn,'Q'))
end

plot([100 130; 100 100]',[100 100; 100 130]','--');
legend('no flow','constant head','well','image well')


% AXI-SYMMETRIC MODEL

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(0);

% set grid
m.setgrid(logspace(-1,7,81),[1;10;10],true);

% set parameters
m.par.constant = [true;false(2,1)];
m.par.kr = [0;10;5];
m.par.cz = [1000;500];

% set stresses
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(2:3,1) = [750;500];

% set solver
m.setsolver(1e-5,20,5);

% run model
m.run;


% SUPERPOSITION MODEL

% 2D raster
[x,y] = meshgrid(0:100);

% interpolation and superposition
s = zeros([size(x),m.grid.nz-1]);
for j = 1:length(xp)
    for i = 1:m.grid.nz-1
        r = sqrt((x-xp(j)).^2+(y-yp(j)).^2);
        tmp = reshape(m.interp1r(i+1,r(:)),size(r));
        s(:,:,i) = s(:,:,i) + q(j) * tmp;
    end
end

% contour plot
figure
for i = 1:2
    subplot(1,2,i)
    contourf(x,y,s(:,:,i))
    colorbar
    axis equal
    xlim([0 100])
    ylim([0 100])
    set(gca,'fontsize',12)
    xlabel('x (m)')
    ylabel('y (m)')
    title(['aquifer ',int2str(i)])
end

% effective velocities
n = 0.3;
[x,y] = deal(x(2:end-1,2:end-1),y(2:end-1,2:end-1));
ip = find(x==xp(1) & y==yp(1));
di = [0 1 -1 98:100 -100:-98];
for i = 1:size(s,3)
    tmp = -m.par.kr(i+1)/n*diff(s(:,:,i),1,2);
    tmp = (tmp(2:end-1,1:end-1) + tmp(2:end-1,2:end))/2;
    tmp(ip+di) = 0;
    qx(:,:,i) = tmp;
    tmp = -m.par.kr(i+1)/n*diff(s(:,:,i),1,1);
    tmp = (tmp(1:end-1,2:end-1) + tmp(2:end,2:end-1))/2;
    tmp(ip+di) = 0;
    qy(:,:,i) = tmp;
end

% quiver plot
figure
for i = 1:2
    subplot(1,2,i)
    quiver(x,y,qx(:,:,i)/10,qy(:,:,i)/10,0)
    axis equal
    xlim([80 100])
    ylim([80 100])
    set(gca,'fontsize',12)
    xlabel('x (m)')
    ylabel('y (m)')
    title(['aquifer ',int2str(i)])
end