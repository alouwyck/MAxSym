% example of slug interference test in multi-layer system
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set grid
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2), rb];
m.setgrid(rb,[2;0.5;0.5;1;0.5;0.5;1],false);

% set parameters
m.par.inactive = false(m.grid.nz,m.grid.nr);
m.par.inactive([1:3 5:7],1) = true;
m.par.kr = [5;ones(6,1)];
m.par.cr = [0.01,zeros(1,m.grid.nr-2)];
m.par.kz = [1;0.1*ones(6,1)];
m.par.cz = [10;zeros(5,1)];
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-5*ones(1,m.grid.nr-1)];
m.par.sy = 0.2;

% set stresses
m.stress.s0 = sparse(m.grid.nz,m.grid.nr);
m.stress.s0(4,1) = 2;

% set solver
m.setsolver(1e-6,20,5);

% run model
m.run;

% drawdown at 1m
s1 = squeeze(m.interp1r(4,1,[]));

% time-drawdown graph
figure
semilogx(m.time.t,[squeeze(m.s(4,1,:)),s1],'-')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 0.1])
legend('slugged well','1.00m')

% plot model set up
figure
x = [0.01 m.grid.rb(end)];
y = m.grid.zb';
set(axes('parent',gcf),'xlim',x,'xscale','log',...
    'ylim',y([end 1]),'box','on','layer','top','fontsize',12)
xlabel('r (m)')
ylabel('z (m)')

% layers and sub-layers
patch([x fliplr(x)],y([1 1 end end]),0.5*ones(1,3),'parent',gca)
rw = m.grid.rb(2);
hl = line([rw,x(end)],[y;y],'color','k','parent',gca);
set(hl(2),'linewidth',5) % clay layer

% well
patch([x(1),rw,rw,x(1)],y([1 1 end end]),'w','parent',gca)
line([rw rw],y([1 end]),'color','k','linewidth',4,'parent',gca) % skin
y = linspace(y(4),y(5),10);
line([x(1),rw],[y;y],'color','k','parent',gca) % screen

% observation well
dx = 0.1;
patch(1+dx*[-1,1,1,-1],[y([end end]),m.grid.zb([1 1])'],'w','parent',gca)
line([1-dx,1+dx],[y;y],'color','k','parent',gca) % screen
