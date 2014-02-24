% example of slug test interpretation
% A.LOUWYCK (2011)


% SYNTHETIC SLUG TEST

% MAxSym model
m = MAxSym.Model;
m.settime(diff(logspace(-6,1,251)));
rb = logspace(log10(0.03),3,401);
rb = [rb(1)^2/rb(2),rb];
m.setgrid(rb,1,true);
m.par.kr = 1;
m.par.ss = [rb(2)^2*pi/m.grid.vol(1),1e-5*ones(1,m.grid.nr-1)];
m.stress.s0 = [1,zeros(1,m.grid.nr-1)];
m.setsolver(1e-5,1);
m.run;

% observations
tobs = (1:1e4)/1440/60;
sobs = squeeze(m.interp1t(1,1,tobs));
sobs = sobs + 0.001*randn(size(sobs));
sobs(sobs<0) = 0;


% INTERPRETATION

% regression parameters
dp = 10^0.1;
delta = 1e-5;
mni = 50;
mu = 10;

% initial run
n = 0;
m.par.kr = 10;
m.par.ss(2:end) = 1e-3;
m.run;
s = squeeze(m.interp1t(1,1,tobs));
eta = sobs(:)-s(:);
ssr = eta' * eta;
dssr = delta + 1;

% echo
fprintf(1,'\nInitial run\n');
fprintf(1,' T = %e\n',m.par.kr)
fprintf(1,' S = %e\n',m.par.ss(end))    
fprintf(1,' SSR = %e\n',ssr);

% iterations
while n < mni && dssr > delta
    
    % sensitivities
    m.par.kr = m.par.kr*dp;
    m.run;
    tmp = squeeze(m.interp1t(1,1,tobs));
    J(:,1) = (tmp(:)-s(:))/log10(dp);
    m.par.kr = m.par.kr/dp;
    
    m.par.ss(2:end) = m.par.ss(2:end)*dp;
    m.run;
    tmp = squeeze(m.interp1t(1,1,tobs));
    J(:,2) = (tmp(:)-s(:))/log10(dp);
    m.par.ss(2:end) = m.par.ss(2:end)/dp;
    
    % condition number
    [~,v] = svd(J);
    v = diag(v);
    k = max(v)/min(v);

    % adjusting parameters
    H = J'*J;
    B = 10.^((H+mu*eye(size(H)))\(J'*eta));
    m.par.kr = m.par.kr*B(1);
    m.par.ss(2:end) = m.par.ss(2:end)*B(2);
    
    % residuals and sum of squares
    m.run;
    s = squeeze(m.interp1t(1,1,tobs));
    eta = sobs(:)-s(:);
    dssr = eta' * eta;
    decrease = dssr < ssr;
    [dssr,ssr] = deal(abs((ssr-dssr)/ssr),dssr);
    
    % iteration index
    n = n + 1;
    
    % echo
    fprintf(1,'\nIteration %d\n',n)
    fprintf(1,' condition number = %.2f\n',k)
    fprintf(1,' T = %e\n',m.par.kr)
    fprintf(1,' S = %e\n',m.par.ss(end))    
    fprintf(1,' SSR = %e\n',ssr);
    fprintf(1,' damping factor = %e\n',mu);
    
    % damping factor
    if decrease
        mu = mu/2;
    else
        mu = mu*2;
    end
    
end

% time-drawdown plot
figure
semilogx(tobs,sobs,'k-',tobs,s,'r-')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')

% sensitivity plot
figure
semilogx(tobs,J)
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('sensitivity')
legend('T','S')
