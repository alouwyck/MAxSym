% example of model for transient well field in multi-aquifer system with recharge
% A.LOUWYCK (2011)


% MAXSYM MODELS

% recharge
N = 5e-4;
R = sqrt(250*20/pi/N);

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,3,81)));

% set grid
rb = logspace(-5,log10(R),51);
rb = [rb(1:end-1) logspace(log10(R),7,41)];
m.setgrid(rb,[20;10;30],false);

% set parameters
m.par.kr = [20;0.1;10];
m.par.kz = [10;0.01;1];
m.par.ss = [1e-4;1e-6;1e-5];
m.par.sy = 0.2;

% set solver
m.setsolver(1e-7,50,5);

% run model with recharge
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(1,2:50) = -N*m.grid.hs(2:50);
m.run;
sr = permute(m.s,[2 1 3]);

% run model with pumping
m.stress.q = zeros(m.grid.nz,m.grid.nr);
m.stress.q(3,1) = 1;
m.run
sw = permute(m.s,[2 1 3]);


% CONTOUR PLOT

% wells
[xp,yp,q] = deal(zeros(20,1),-95:10:95,250*ones(20,1));

% 2D raster
[x,y] = meshgrid(-500:500);

% interpolation and superposition: pumping wells
ilay = 3;
it = [41 81];
s = zeros([size(x),length(it)]);
for j = 1:length(xp)
    r = sqrt((x(:)-xp(j)).^2+(y(:)-yp(j)).^2);
    tmp = interp1(log10(m.grid.r),squeeze(sw(:,ilay,it)),log10(r));
    b = r < m.grid.r(1);
    if any(b)
        tmp(b,:) = squeeze(repmat(sw(1,ilay,it),[sum(b) 1 1]));
    end
    s = s + q(j) * reshape(tmp,size(s));
end

% interpolation and superposition: recharge
r = sqrt(x(:).^2+y(:).^2);
tmp = interp1(log10(m.grid.r),squeeze(sr(:,ilay,it)),log10(r));
b = r < m.grid.r(1);
tmp(b,:) = squeeze(repmat(sr(1,ilay,it),[sum(b) 1 1]));
s = s + reshape(tmp,size(s));

% plot
figure
for i = 1:length(it)
    subplot(1,length(it),i)
    contourf(x,y,s(:,:,i))
    colorbar
    caxis([-10 0])
    axis equal
    xlim([-500 500])
    ylim([-500 500])
    set(gca,'fontsize',12)
    xlabel('x (m)')
    ylabel('y (m)')
    title([num2str(m.time.t(it(i)),'%.0e'),' days of pumping'])
end


% TIME-DRAWDOWN PLOT

% point
[x,y] = deal(0);

% interpolation and superposition: pumping wells
s = zeros(m.grid.nz,length(m.time.t));
for j = 1:length(xp)
    r = sqrt((x-xp(j)).^2+(y-yp(j)).^2);
    if r > m.grid.r(1)
        tmp = interp1(log10(m.grid.r),sw,log10(r));
    else
        tmp = sw(1,:,:);
    end
    s = s + q(j) * squeeze(tmp);
end

% interpolation and superposition: recharge
r = sqrt(x.^2+y.^2);
if r > m.grid.r(1)
    tmp = interp1(log10(m.grid.r),sr,log10(r));
else
    tmp = sr(1,:,:);
end
s = s + squeeze(tmp);

% plot
figure
semilogx(m.time.t,s')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
legend('upper aquifer','aquitard','lower aquifer')

