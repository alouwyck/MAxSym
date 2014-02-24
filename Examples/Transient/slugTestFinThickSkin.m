% example of slug test in well with finite-thickness skin
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-6,1,251)));

% set grid
rb = logspace(log10(0.03),log10(0.05),20);
rb = [rb(1)^2/rb(2), rb(1:end-1), logspace(log10(0.05),3,381)];
m.setgrid(rb,1,true);

% set parameters
m.par.kr = ones(1,m.grid.nr);
m.par.ss = [m.grid.rb(2)^2*pi/m.grid.vol(1), 1e-5*ones(1,m.grid.nr-1)];

% set stresses
m.stress.s0 = [1 ,zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% loop through different kr values for skin
kskin = logspace(-2,2,5);
for n = 1:length(kskin)

    % set kr
    m.par.kr(1:20) = kskin(n);
    
    % run model
    m.run;
    
    % get well drawdown
    s(:,n) = squeeze(m.s(1,1,:));
    
end

% time-drawdown graph
figure
semilogx(m.time.t,s)
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
xlim([1e-6 1])
legend(strcat(num2str(kskin','%g'),'m²/d'))