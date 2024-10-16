clc;  beep off; close all; clear all; curr_dir = cd;

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

p.vstar = 0.1;
p.l = 0.02688; % always fixed
p.wx = 0.05; % wind strength

% -------------------------------------------------------------------------
% get the fixed points
% -------------------------------------------------------------------------

cd('../code');

u = linspace(1e-3,3,1000);
Y = p.vstar*exp(-1./u)-p.l*u;
% plot(u,Y,'LineWidth',2);

index = find(Y(1:end-1).*Y(2:end) < 0);

fun = @(u) p.vstar*exp(-1./u)-p.l*u;

fp1 =  find_zero(u(index(1)),u(index(1)+1),fun);
fp2 = find_zero(u(index(2)),u(index(2)+1),fun);

s.UL = [0;0];
s.UR = [fp2;0];

% -------------------------------------------------------------------------
% solve the profile
% -------------------------------------------------------------------------

N = 30;
alpha_vals = linspace(1e-3,1-1e-3,N);
eig_vals = zeros(size(alpha_vals));
c_vals = zeros(size(alpha_vals));

for j=1:length(alpha_vals)

    p.alpha = alpha_vals(j);
     
    [p,s] = get_profile(p,s);
    
    fun_acc = @(omega)get_acc_angle(omega,s,p);

    eig_vals(j) = find_zero(-0.01,0.01,fun_acc);

    c_vals(j) = p.c;

end

%plot_profile(s);

plot(c_vals,eig_vals,'-k','LineWidth',2);
h = xlabel('c');
set(h,'FontSize',18);
h = ylabel('eigenvalue');
set(h,'FontSize',18);

cd(curr_dir);