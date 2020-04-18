beep off; clc; clear all; close all; curr_dir = cd;
addpath(strcat(pwd,'/../bin/'));

%
% parameters
%

p.a = 4; % a > 0
p.gamma = 2; % gamma >= 1

%
% dependent variables
%

p.mu = (p.a-sqrt(p.gamma^2-1))/(p.a+sqrt(p.gamma^2-1));
p.nu = 1/(p.a+sqrt(p.gamma^2-1));

%
% structure variables
%

s.I=300; 
s.R=s.I;
s.L=-s.I;

%
% solve for the eigenfunction
%

L = s.I;
eigfun_ode = @(x,y,lambda)([A(x,lambda,s,p)*y(1:4); ...
    A(x+L,lambda,s,p)*y(5:8)]);

fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p));

guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];

lambda0 = 0.0586 + 1.3043i;

guess = @(x)[guess_fun(x);guess_fun(x+L)];

solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',1e-10,'AbsTol',1e-10,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters;

%% plot eigenfunction
x = linspace(-L,0,500);
y = deval(eigfun_sol,x);
hold on;
plot(x,real(y(1:4,:)),'-g','LineWidth',2)
plot(x+L,real(y(5:8,:)),'-b','LineWidth',2)


% degree of Chebyshev interpolant of the solution of the homological
% equation
N = 4001;

% truncated domain endpoints
a = -L; 
b = L;

% Transformation to get Chebyshev coefficients (to speed up the
% computation, we use a Chebyshev interpolant of the solutions to the
% homological equations).

% Chebyshev nodes
theta = ((0:1:N-1)+0.5)*pi/N;
x0 = cos(theta); % in [-1,1]
x = 0.5*(a+b)+0.5*(a-b)*x0; % in [a,b]

Id2 = (2/N)*speye(N);
Id2(1,1) = Id2(1,1)/2;
Tcf = Id2*cos(theta.'*(0:1:N-1)).';

cf = Tcf*eval_eigfun(eigfun_sol,x,s.I).';

data.lam = lam;
data.cf = cf;
data.N = N;
data.a = a;
data.b = b;
data.sol = eigfun_sol;
data.p = p;

cd('../data');
save('eig_fun_Schrodinger_a_4_gamma_2','data');
cd(curr_dir);


