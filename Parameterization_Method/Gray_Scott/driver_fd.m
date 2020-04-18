clc; clear all; close all; beep off; curr_dir = cd;

% load the solutions to the homological equations
addpath(strcat(pwd,'/../bin/'));
addpath(strcat(pwd,'/../bin/dmsuite'));
cd('../data');
ld = load('P_Gray_Scott');
clc;
data = ld.data;
p = data.p;
P = data.P; % structure containing solutions to the homological equations
L = data.L; % Numerical infinity for homological equation solutions
% num_homological_equations = 5; % number of homological equations solved
lam = data.lam; % eigenvalue for Nagumo's equation
cd(curr_dir);

desired_error_tol = 1e-8;

figure;
hold on;
sigma0 = 0.1;
sigma1 = (desired_error_tol/data.TOL)^(1/(length(P)-1))
T = (log(sigma1)-log(sigma0))/lam

return
xgrid = linspace(-L,L,1001);
u0 = evaluate_homological_equations(P,sigma0,xgrid);
u1 = evaluate_homological_equations(P,sigma1,xgrid);

hold on;
plot(xgrid,u0,'-k','LineWidth',2);
plot(xgrid,u1,'-g','LineWidth',2);

drawnow;
% return

%% User defined FD parameters

for H = [2^-2,2^-3,2^-4,2^-5,2^-6,2^-7,2^-8] % distance between spatial nodes

    constant = 0.25; % (Delta t)/(Delta x)^2 = constant
    L_fd = data.L; % To verify the conjugacy, we use time evolution restricted to the interval [-L_fd,L_fd].
    stats = 'off'; % print and plot information about the time evolution

    t0 = tic;
    [U_n,u0,u1] = fd_Gray_Scott(p,P,T,H,constant,sigma0,L_fd,data.lam,stats);
    t1 = toc(t0);

    % The difference at the final time between the parameterization
    % solution and the finite difference solution.
    
    diff = max(max(abs(u1([1,3],:)-U_n)));
    fprintf('Delta x = %4.4g, run time = %4.4g, difference = %4.4g\n\n',H,t1,diff);
    
    
end


