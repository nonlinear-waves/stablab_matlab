function plot_profile_local_model(p,s)

x = linspace(s.L,s.R,200);
y = zeros(length(x),2);
for j = 1:length(x)
	temp = soln(x(j),s);
	y(j,1) = temp(1);
	y(j,2) = temp(2);
end

figure;
hold on;


% wave speed
spd = p.spd;

% values at - infinity
v_neg = p.v_neg;
tau_neg = p.tau_neg;
S_neg = p.S_neg;

% parameters
m1 = p.m1;
none = p.none;
mu = p.mu;
kappa = p.kappa;
r = zeros(length(x),1);
w = zeros(length(x),1);

for j = 1:length(x)
tau = y(j,1);
S = y(j,2);
v = v_neg - spd * (tau - tau_neg);
r(j,1) = v;
T = exp(S) / tau + 0.1e1;
w(j,1) = T;
end

plot(x,r,x,y,x,w,'LineWidth',2);

h = gca;
set(h,'FontSize',22);
set(h,'LineWidth',2);

h = xlabel('x');
set(h,'FontSize',22);
h = ylabel('profiles');
set(h,'FontSize',22);

h = legend('v','tau','S','T','Location','Best');
set(h,'FontSize',22);
