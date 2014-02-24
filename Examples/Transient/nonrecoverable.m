% example of recovery test with non-recoverable compression 
% A.LOUWYCK (2011)


% PUMPING test

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(diff(logspace(-5,0,51)));

% set grid
m.setgrid(logspace(-1,7,81),1,true);

% set parameters
m.par.kr = 10;
m.par.ss = 1e-3;

% set stresses
m.stress.q = [100, zeros(1,m.grid.nr-1)];

% set solver
m.setsolver(1e-5,1);

% run model
m.run;

% get drawdown matrix
s = squeeze(m.s);


% RECOVERY phase

% modify stresses
m.stress.q = [];
m.stress.s0 = s(:,end)';

% loop through different Ss values
ss = 10.^(-5:-3);
for n = 1:length(ss)

    % set ss
    m.par.ss = ss(n);

    % run model
    m.run
    
    % get drawdown
    sr(:,:,n) = squeeze(m.s);

end

% time-drawdown graph
figure
ir = 11;
t = m.time.t(:);
s = cat(2,repmat(s,[1 1 length(ss)]),sr(:,2:end,:));
plot([t; t(end)+t(2:end)],squeeze(s(ir,:,:)))
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
title(strcat(num2str(m.grid.r(ir)),'m'));
legend(num2str(ss','%.0e'))
