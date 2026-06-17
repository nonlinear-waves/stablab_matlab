clc; close all; clear all; beep off;

% test that mu = 0 does not correspond to an eigenvalue

% parameters
p.rho_m = 0.1;
p.rho_p = 1; % bigger than p.rho_m
p.h = 1.5;
p.L = 1;

% don't change
p.g = 9.8;
p.a = log(p.rho_p/p.rho_m)/p.h;


ode = @(x,y)[0,1; -p.a^2/4,-p.a]*y;
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol = ode15s(ode,[-p.h,0],[0;1],options);

temp = deval(sol,0);

constant = -p.g*temp(2)/temp(1)

k2 = -p.a*temp(2)/temp(1)-p.a^2/4

term = (2/p.h)*(2*p.h-p.a^2)-1/(4*p.h^2)