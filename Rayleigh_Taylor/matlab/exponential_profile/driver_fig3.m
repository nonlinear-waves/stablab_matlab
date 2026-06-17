clc; close all; clear all; beep off;

disp('Plot of the top 5 modes for k = 200');

% parameters
p.rho_m = 0.1;
p.rho_p = 1; % bigger than p.rho_m
p.h = 1.5;
p.L = 1;

% don't change
p.g = 9.8;
p.a = log(p.rho_p/p.rho_m)/p.h;
profile_ratio = @(x)p.a;
s.options = odeset('RelTol',1e-8,'AbsTol',1e-8);

Atwood_number =  p.rho_p -p.rho_m/(p.rho_p + p.rho_m);
fprintf('\nAtwood number = %4.4g\n',Atwood_number);

bound = sqrt(p.a*p.g);

% compute evans function

k = 200;

rts = zeros(1,5);
index = 1:1:5;
for j = 1:length(rts)
    rts(j) = get_eigenvalues(index(j),k,p.a,p.h);
end

figure;
hold on;
plot([1,5],[bound,bound],'--b','LineWidth',2);
plot(rts,'.-k','MarkerSize',18,'LineWidth',2);
axis([1,5,3.86,1.001*max(rts)]);
h = xlabel('Mode index j');
set(h,'FontSize',18);
h = ylabel('\lambda');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);