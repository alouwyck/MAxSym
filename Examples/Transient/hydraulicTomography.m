% example of hydraulic tomography test in multi-layer system
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
dtpump = diff(logspace(-6,log10(20/1440),101));
dtrecov = diff(logspace(-6,log10(100/1440),501));
m.settime(repmat([dtpump,dtrecov],1,3),repmat([100 500],1,3));

% set grid
rb = logspace(log10(0.05),log10(0.07),10);
rb = [rb(1)^2/rb(2), rb];
rb = [rb(1:end-1),logspace(log10(0.07),3,241)];
m.setgrid(rb,ones(3,1),true);

% set parameters
m.par.kr = repmat([10;2;5],1,m.grid.nr);
m.par.kr(:,1:10) = 50;
m.par.kz = repmat([5;0.5;2],1,m.grid.nr);
m.par.kz(:,1:10) = 50;
m.par.cz = sparse(m.grid.nz-1,m.grid.nr);
m.par.cz(:,1) = 1e5;
m.par.ss = repmat([1e-4;1e-5;5e-5],1,m.grid.nr);
m.par.ss(:,2:10) = 1e-6;
m.par.ss(:,1) = m.grid.rb(2)^2*pi./m.grid.vol(:,1);

% set stresses
for n = 1:2:m.time.nper
    m.stress(n).q = sparse(m.grid.nz,m.grid.nr);
    m.stress(n).q((n+1)/2,1) = 10;
end

% set solver
m.setsolver(1e-6,20,5);

% run model
m.run;

% time-drawdown graphs
figure
ir = [1 89];
for n = 1:length(ir)
    subplot(1,length(ir),n)
    plot(m.time.t*24,squeeze(m.s(:,ir(n),:))')
    set(gca,'fontsize',12)
    xlabel('time (h)')
    ylabel('drawdown (m)')
    if n == 1
        title('extraction well')
    else
        title(['obs.well at ',num2str(m.grid.r(ir(n)),'%.2f'),'m'])
    end
    legend(strcat('layer',num2str((1:m.grid.nz)')),'location','southeast')
end
