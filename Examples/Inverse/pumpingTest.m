% example of pumping test interpretation
% A.LOUWYCK (2011)


% SYNTHETIC PUMPING TEST

% MAxSym model
m = MAxSym.Model;
m.settime(diff(logspace(-5,0.1,101)));
m.setgrid(logspace(-1,7,81),1,true);
m.par.kr = 100;
m.par.ss = 1e-3;
m.stress.q = [100, zeros(1,m.grid.nr-1)];
m.setsolver(1e-5,1);
m.run;

% observations
[tobs,robs] = deal((1:1440)/1440,[1 5 10]);
sobs = squeeze(m.interp2(1,robs,tobs))';
sobs = sobs + 0.001*randn(size(sobs)); 


% INTERPRETATION

% regression parameters
dp = 10^0.1;
delta = 1e-5;
mni = 50;

% initial run
n = 0;
m.par.kr = 10;
m.par.ss = 1e-4;
m.run;
s = squeeze(m.interp2(1,robs,tobs))';
eta = sobs(:)-s(:);
ssr = eta' * eta;
dssr = delta + 1;

% echo
fprintf(1,'\nInitial run\n');
fprintf(1,' T = %e\n',m.par.kr)
fprintf(1,' S = %e\n',m.par.ss)    
fprintf(1,' SSR = %e\n',ssr);

% iterations
while n < mni && dssr > delta
    
    % sensitivities
    m.par.kr = m.par.kr*dp;
    m.run;
    tmp = squeeze(m.interp2(1,robs,tobs))';
    J(:,1) = log10(tmp(:)./s(:))/log10(dp);
    J(:,1) = (tmp(:)-s(:))/log10(dp);
    m.par.kr = m.par.kr/dp;
    
    m.par.ss = m.par.ss*dp;
    m.run;
    tmp = squeeze(m.interp2(1,robs,tobs))';
    J(:,2) = log10(tmp(:)./s(:))/log10(dp);
    J(:,2) = (tmp(:)-s(:))/log10(dp);
    m.par.ss = m.par.ss/dp;
    
    % condition number
    [~,v] = svd(J);
    v = diag(v);
    k = max(v)/min(v);

    % adjusting parameters
    B = 10.^((J'*J)\(J'*eta));
    m.par.kr = m.par.kr*B(1);
    m.par.ss = m.par.ss*B(2);
    
    % residuals and sum of squares
    m.run;
    s = squeeze(m.interp2(1,robs,tobs))';
    eta = sobs(:)-s(:);
    dssr = eta' * eta;
    [dssr,ssr] = deal(abs((ssr-dssr)/ssr),dssr);
    
    % iteration index
    n = n + 1;
    
    % echo
    fprintf(1,'\nIteration %d\n',n)
    fprintf(1,' condition number = %.2f\n',k)
    fprintf(1,' T = %e\n',m.par.kr)
    fprintf(1,' S = %e\n',m.par.ss)    
    fprintf(1,' SSR = %e\n',ssr);    
    
end

% plot
figure
plot(tobs',sobs,'k-',tobs',s,'r-')
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')

