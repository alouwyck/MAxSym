% example of model for extraction near infiltration pond
% A.LOUWYCK (2011)


% INPUT

% wells and infiltration pond
[x1,y1,Q1,d1] = deal(-15,0,500,[1 2]/3);
[x2,y2,Q2,d2,f2] = deal(15,0,1000,[1 1]/2,[1.5 0.5]);
N = 0.5;
R = 10;

% aquifer
D = [20;30];
kr = [10;20];
cz = 1000;
ss = [1e-4;1e-5];
sy = 0.2;

% number of days
days = (3*365)+1;

% map and discharge plot
figure
subplot(121)
plot(x1,y1,'ko',x2,y2,'ko')
text(x1,y1,'  Q1')
text(x2,y2,'  Q2')
hold on
a = linspace(0,2*pi,100);
plot(R*cos(a),R*sin(a),'b-')
axis equal
xlim([-20 20])
ylim([-20 20])
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('y (m)')

subplot(122)
h = plot([0 d1([1 1]) 1]*24,[Q1 Q1 0 0],[0 d2([1 1]) 1]*24,Q2*f2([1 1 2 2]));
set(h,'linewidth',3)
xlim([0 24])
ylim([-500 2000])
grid minor
set(gca,'fontsize',12)
xlabel('time (h)')
ylabel('Q (m³/d)')
legend('Q1','Q2')


% MAXSYM MODELS

% recharge model
mr = MAxSym.Model;
mr.settime(diff(logspace(-5,log10(days),81)));
rb = logspace(-1,log10(R),21);
rb = [rb(1:end-1),logspace(log10(R),7,61)];
mr.setgrid(rb,D,false);
mr.par.kr = kr;
mr.par.cz = cz;
mr.par.ss = ss;
mr.par.sy = sy;
mr.stress.q = sparse(mr.grid.nz,mr.grid.nr);
mr.stress.q(1,1:20) = -N*mr.grid.hs(1:20);
mr.setsolver(1e-7,50,5);
mr.run;

% well 1 model
m1 = MAxSym.Model;
t = 0:8:(24*days);
t(1) = 1e-5;
t(3:3:end) = [];
dt = [];
for i = 1:length(t)-1
    dt = [dt, diff(logspace(log10(t(i)),log10(t(i+1)),51))];
end
dt = dt/24;
m1.settime(dt,repmat([50 50],1,days));
m1.setgrid(logspace(-1,7,81),D,false);
m1.par.kr = kr;
m1.par.cz = cz;
m1.par.ss = ss;
m1.par.sy = sy;
for k = 1:2:m1.time.nper
    m1.stress(k).q = sparse(m1.grid.nz,m1.grid.nr);
    m1.stress(k).q(2,1) = Q1;
end
m1.setsolver(1e-7,50,5);
m1.run;

% well 2 model
m2 = MAxSym.Model;
t = 0:12:(24*days);
t(1) = 1e-5;
dt = [];
for i = 1:length(t)-1
    dt = [dt, diff(logspace(log10(t(i)),log10(t(i+1)),51))];
end
dt = dt/24;
m2.settime(dt,repmat([50 50],1,days));
m2.setgrid(logspace(-1,7,81),D,false);
m2.par.kr = kr;
m2.par.cz = cz;
m2.par.ss = ss;
m2.par.sy = sy;
for k = 1:2:m2.time.nper
    m2.stress(k).q = sparse(m2.grid.nz,m2.grid.nr);
    m2.stress(k).q(2,1) = f2(1)*Q2;
end
for k = 2:2:m2.time.nper
    m2.stress(k).q = sparse(m2.grid.nz,m2.grid.nr);
    m2.stress(k).q(2,1) = f2(2)*Q2;
end
m2.setsolver(1e-7,50,5);
m2.run;


% TIME-DRAWDOWN GRAPHS

% input and initialization
[x,y] = deal([x1;x2;0],[y1;y2;0]);
[xb,yb] = deal(x,y);
t = (1:(24*(days-1)))/24;
s = zeros(2,3,length(t));
m = {m1; m2; mr};

% interpolation and superposition
for n = 1:length(m)
    r = sqrt((x-xb(n)).^2+(y-yb(n)).^2);
    s = s + m{n}.interp2([],r(:),t(:));
end
s = reshape(permute(s,[3 2 1]),[],6);

% plots
figure
b = 8:24:length(t);
plot(t(b),s(b,:))
hold on
b = 24:24:length(t);
plot(t(b),s(b,:))
xlim(t([1 end]))
set(gca,'fontsize',12)
xlabel('time (d)')
ylabel('drawdown (m)')
str = sprintf('layer %1d: (%d,%d);',[[ones(3,1);2*ones(3,1)],repmat([x(:) y(:)],2,1)]');
str = regexp(str,';','split');
legend(str(1:end-1),'location','northeastoutside');

figure
b = t >= days-2;
plot((t(b)-days+2)*24,s(b,:),'-x')
xlim([0 24])
set(gca,'fontsize',12)
xlabel('time (h)')
ylabel('drawdown (m)')
grid minor
legend(str(1:end-1),'location','northeastoutside');


% DISTANCE-DRAWDOWN PLOT

% input and initialization
[x,y] = meshgrid(-1000:1000,0);
t = (days-1) - [d1(2) 0];
s = zeros(2,numel(x),2);

% interpolation and superposition
for n = 1:length(m)
    r = sqrt((x-xb(n)).^2+(y-yb(n)).^2);
    b = r > m{n}.grid.r(1);
    s = s + m{n}.interp2([],r(:),t(:));
end
s = reshape(permute(s,[2 3 1]),[size(x) 4]);

% plot
figure
plot(x,squeeze(s))
set(gca,'fontsize',12)
xlabel('x (m)')
ylabel('drawdown (m)')
str = sprintf('layer %1d: %d hours;',[[1 1 2 2];repmat([d1(1) 1]*24,1,2)]);
str = regexp(str,';','split');
legend(str(1:end-1),'location','southeast');

