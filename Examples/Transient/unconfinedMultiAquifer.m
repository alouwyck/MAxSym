% example of pumping test in unconfined multi-aquifer system
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,2,71)));

% set grid
rb = logspace(log10(0.25),log10(0.3),11);
rb = [rb(1)^2/rb(2), rb];
rb = [rb(1:end-1),logspace(log10(0.3),6,71)];
m.setgrid(rb,[2;3;1;1;1;2;1;1;3;4],false);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive(1:3,1:11) = true;
m.par.inactive([4:5,7:8],1) = true;
m.par.kr = repmat([10;10;0.5;5;5;5;5;5;0.01;1],1,m.grid.nr);
m.par.kr(4:8,1:11) = 50;
m.par.kz = repmat([1;1;0.1;1;1;1;1;1;0.0005;0.2],1,m.grid.nr);
m.par.kz(4:8,1:11) = 50;
m.par.ss = repmat(1e-5*[10;10;10;1;1;1;1;1;0.1;1],1,m.grid.nr);
m.par.ss(6,1) = m.grid.rb(2)^2*pi/m.grid.vol(1);
m.par.sy = 0.2;

% set stresses
m.stress.q = sparse(m.grid.nz,m.grid.nr);
m.stress.q(6,1) = 250;

% set solver
m.setsolver(1e-6,50,5);

% run model
m.run;

% time-drawdown graph
figure
ilay = [1 6 10];
ir = 19;
semilogx(m.time.t,[squeeze(m.s(6,1,:)),squeeze(m.s(ilay,ir,:))'])
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
str = {'pumping well'};
for n = ilay
    str{end+1} = sprintf('layer %d, r = %.2fm',n,m.grid.r(ir));
end
legend(str)

% plot model set up
figure
x = [0.1 m.grid.rb(end)];
y = m.grid.zb';
set(axes('parent',gcf),'xlim',x,'xscale','log',...
    'ylim',y([end 1]),'box','on','layer','top','fontsize',12)
x = [x fliplr(x)];
xlabel('r (m)')
ylabel('z (m)')

% layers
nl = cumsum([1 2 1 5 1 1]);
c = {'y','g','y','g','y'};
for n = 1:5
    patch(x,y(nl([n,n,n+1,n+1])),c{n},'parent',gca)
end

% well
rw = m.grid.rb(2);
rs = m.grid.rb(12);
patch([x(1),rw,rw,x(1)],[y([1 1]),y(nl([4 4]))],'b','parent',gca)

% clay seal and skin
patch([rw,rs,rs,rw],[y([1 1]),y(nl([3 3]))],'k','parent',gca)
patch([rw,rs,rs,rw],y(nl([3 3 4 4])),'r','parent',gca)

% sub-layers
line([rw,m.grid.rb(end)],[y;y],'color','k','parent',gca)

% screen
y = linspace(y(6),y(7),10);
line([x(1),rw],[y;y],'color','k','parent',gca)
