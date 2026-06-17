clc; close all; clear all; beep off;

disp('Plot of the largest 5 eigenvalues for k = 680');

% parameters
p.rho_m = 0.1;
p.rho_p = 0.2; % bigger than p.rho_m
p.h = 2;
p.L = 1;

% fixed
p.g = 9.8;
p.rho_0_der = (p.rho_p-p.rho_m)/p.h;
rho_0 = p.rho_m +p.rho_0_der*(linspace(-p.h,0,1000)+p.h);
p.L0 = 1/max(p.rho_0_der./rho_0);

choices = 1:1:5; % 1 for top eigenvalue, 2 for second top eigenvalue
k_max = 680;

left = 1.5;
right = 2.193;


p.k = k_max/p.L;

rts = get_eigenvalues(p,choices,left,right,1000);

figure;
hold on;
plot([1,5],[sqrt(p.g/p.L0),sqrt(p.g/p.L0)],'--b','LineWidth',2);
plot(rts,'.-k','MarkerSize',18,'LineWidth',2);
axis([1,5,0,1.1*max(rts)]);
h = xlabel('Mode index j');
set(h,'FontSize',18);
h = ylabel('\lambda');
set(h,'FontSize',18);
h = gca;
set(h,'FontSize',18);




