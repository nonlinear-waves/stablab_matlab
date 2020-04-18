clc; clear all; close all; beep off; curr_dir = cd;

% load the solutions to the homological equations
addpath(strcat(pwd,'/../bin/'));
addpath(strcat(pwd,'/../bin/dmsuite'));

cd('../data');
ld = load('P_Schrodinger');
clc;
data = ld.data; 
p = data.p; % system parameters
P = data.P; % structure containing solutions to the homological equations
L = data.L; % Numerical infinity for homological equation solutions
num_homological_equations = length(P)-1; % number of homological equations solved
lam1 = data.lam;
cd(curr_dir);

%% FD set up.

theta1 = 1.5;
theta2 = 1.5;


T = 1; % time period for which we do time evolution

constant = 0.25; % (Delta t)/(Delta x)^2 = constant
L_fd = 30; % To verify the conjugacy, we use time evolution restricted to the interval [-L_fd,L_fd].
stats = 'off'; % print and plot information about the time evolution

for H = [2^-2,2^-3,2^-4,2^-5,2^-6,2^-7,2^-8] % distance between spatial nodes


    t0 = tic;
    [U_n,u0,u1] = fd_Schrodinger(p,P,T,H,L_fd,constant,theta1,theta2,lam1,stats);   
    t1 = toc(t0);

    % The difference at the final time between the parameterization
    % solution and the finite difference solution.
    
    diff = max(max(abs(u1-U_n)));
    fprintf('Delta x = %4.4g, run time = %4.4g, difference = %4.4g\n\n',H,t1,diff);
    
    
end




