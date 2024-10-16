clc;  beep off; close all;  clear all; curr_dir = cd;

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

p.vstar = 0.1;
p.l = 0.02688; % always fixed
speed_guess = -0.0854708;

wx = 0.05; % wind strength

N = 30;
alpha_vals = linspace(0,1,N+2);
alpha_vals = alpha_vals(2:end-1);

% -------------------------------------------------------------------------
% get the fixed points
% -------------------------------------------------------------------------

cd('../code');

p.wx = 0;

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
% solve the profile equation for c_star
% -------------------------------------------------------------------------

p.alpha = 0.5;
p.wx = wx;

s.I = 150;
s.scale = 0.1;

guess_fun = @(x)guess(x,s);
s = solve_profile(s,p,speed_guess,guess_fun);
p.cstar = s.sol.parameters;

% plot_profile(s);
% return

dom = linspace(s.L,s.R,4000);
for j = 1:length(alpha_vals)

    p.alpha = alpha_vals(j);
        
    guess_fun = @(x)[soln(x,s);soln(-x,s)];
    
    s = solve_profile(s,p,speed_guess,guess_fun);
    p.c = s.sol.parameters;

    y = zeros(2,length(dom));
    for k = 1:length(dom)
        y(:,k) = soln(dom(k),s);
    end

    mat = [p.wx,p.vstar,p.l,p.wx,p.alpha,p.cstar,p.c];

    name = ['id',num2str(j)];

    cd('../profile_data');
    writematrix(mat,[name,'_info.csv']);
    writematrix(y,[name,'_y.csv'])
    cd('../code');

end
cd(curr_dir);


% hold on;
% dom = linspace(s.L,s.R,4000);
% y = zeros(2,length(dom));
% for j = 1:length(dom)
%     y(:,j) = soln(dom(j),s);
% end
% plot(dom,y,'-r','LineWidth',2);
% h = xlabel('x');
% set(h,'FontSize',18);
% h = ylabel("u,u'");
% set(h,'FontSize',18);
% h = gca;
% set(h,'FontSize',18);


