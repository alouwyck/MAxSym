% example of 2D cross-sectional model to simulate 2 permanent wells in a leaky aquifer
% A.LOUWYCK (2011)

% create Model object
m = MAxSym.Model;

% set time steps
m.settime(0);

% set grid
m.setgrid(logspace(-1,7,81),0.1*ones(111,1),true);

% set parameters
m.par.constant = [true; false(110,1)];
m.par.kr = [0.01*ones(51,1); 10*ones(60,1)];
m.par.kz = m.par.kr/10;

% set solver
m.setsolver(1e-5,100,5);

% run models
su = [];
for n = [70 90]
    
    % discharge
    m.stress.q = zeros(m.grid.nz,m.grid.nr);
    m.stress.q(n+(1:10),1) = 0.1;

    % run
    m.run;
    su = cat(3,su,m.s);

end

% cross section
x = -50:0.1:50;
z = [m.grid.z(:);m.grid.zb(end)];

% wells
[x1,q1] = deal(-2.5,100);
[x2,q2] = deal(2.5,250);
r1 = abs(x-x1);
r2 = abs(x-x2);

% interpolation
for n = 2:m.grid.nz
    s1(n,:) = interp1(log10(m.grid.r),su(n,:,1),log10(r1));
    s1(n,r1<m.grid.r(1)) = su(n,1,1);
    s2(n,:) = interp1(log10(m.grid.r),su(n,:,2),log10(r2));
    s2(n,r2<m.grid.r(1)) = su(n,1,2);
end

% superposition
s = q1*s1 + q2*s2;
s = [s; s(end,:)];

% contour plot
figure
[cs,h] = contourf(x,z,log10(-s),-1.4:0.1:1.3);
clabel(cs,h,-0.2:0.2:0.6)
colorbar
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('z (m)')

% specific discharges
[x,z] = deal(x(2:end-1),z(2:end-1));
qx = -repmat(m.par.kr([1:end,end]),1,size(s,2)-1) .* diff(s,1,2)/0.1;
cv = m.grid.D ./ m.par.kz / 2;
cv = [cv(1:end-1) + cv(2:end); cv(end)];
qz = diff(s,1,1) ./ repmat(cv,1,size(s,2));
qx = (qx(2:end-1,1:end-1) + qx(2:end-1,2:end))/2;
qz = (qz(1:end-1,2:end-1) + qz(2:end,2:end-1))/2;

% specific discharge plot 
figure

subplot(1,2,1)
iz = z<=6;
ix = -2:1:2;
plot(qx(iz,ismember(x,ix)),z(iz))
legend(strcat(num2str(ix'),'m'))
set(gca,'fontsize',12)
xlabel('q_x (m/d)')
ylabel('z (m)')

subplot(1,2,2)
plot(qz(iz,ismember(x,ix)),z(iz))
legend(strcat(num2str(ix'),'m'))
set(gca,'fontsize',12)
xlabel('q_z (m/d)')
ylabel('z (m)')