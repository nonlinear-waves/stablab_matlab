
%% Nagumo - reaction diffusion equation ( U_t = U_{xx} - U + U^3 )

%{
Numericaly solves for the unstable manifold solution to Nagumo's equation
 associated with the unstable traveling wave solution, $\sqrt{2}\sech(x)$. 
%}

addpath(strcat(pwd,'/../bin/'));
addpath(strcat(pwd,'/../bin/dmsuite'));
clc; clear all; close all; beep off; curr_dir = cd;

%% User defined parameters

% Spatial infinity (chosen so the $|p_0(L)| < 1e-10$, where $p_0$ is the profile.
L = 25; 
% Scale of the eigenfunction. Chosen so that the last
% homological equation we solve for is on the order of 1e-15 or smaller.
scl = 0.33; 
% The number of homological equations we solve for.
num_homological_equations = 30; 

%% Set up variables

lam = 3; % eigenvalue of linearized PDE
n_vals = 2:1:num_homological_equations;  % homological equations indices.

%% Find the homological equations, P{n}.fun

% Unstable wave solution 
fprintf('n=0, Loading function\n\n');
P{1}.fun = @(x)([sqrt(2)*sech(x); ...
                sqrt(2)*(-tanh(x).*sech(x))]);

% Eigenfunction
fprintf('n=1, Loading function\n\n');
P{2}.fun = @(x) scl*eig_fun(x);

%
% Solve for the solutions to the homological equations
%

options = bvpset('RelTol', 1e-12, 'AbsTol', 1e-12,'Nmax', 20000);
cnt = 0;
for n = n_vals
    
    fprintf(strcat('n=', num2str(n), ', Generating function\n\n'));
    
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam));
    
    bc_handle = @(ya,yb)(bc_fun(ya,yb,n)); 
    solinit = bvpinit([-L,L],P{n}.fun);
    sol = bvp5c(ode_handle,bc_handle,solinit);
            
    P{n+1}.mu_L = sqrt(3*n+1);
    P{n+1}.mu_R = -sqrt(3*n+1);
    P{n+1}.sol = sol;
    
    % degree of polynomial
    N = 240;
    
    a = -L;
    b = L;

    % Chebyshev nodes
    theta = ((0:1:N-1)+0.5)*pi/N;
    x0 = cos(theta); % in [-1,1]
    x = 0.5*(a+b)+0.5*(a-b)*x0; % in [a,b]

    % Transformation to get Chebyshev coefficients
    Id2 = (2/N)*speye(N);
    Id2(1,1) = Id2(1,1)/2;
    Tcf = Id2*cos(theta.'*(0:1:N-1)).';
    cf = Tcf*deval(sol,x).';
    
    P{n+1}.fun = @(x)(interp_cheby_fun(cf,x,n,N,a,b));
    
    nrm = max(norm(deval(sol,x)))

end

%% Plot the homological equations
fprintf('Preparing graphs\n\n');
x = linspace(-L,L,1000);

hold on;
for j = 1:num_homological_equations
    y = P{j}.fun(x);
    plot(x,y,'LineWidth',2)
    nrm = max(max(abs(y)));
end






%% User defined FD parameters
T = 0.15;
H = 0.025;
constant = 0.25;
sigma0 = 1;
L_fd = 10;

K = H^2*constant;
sigma1 = sigma0*exp(lam*T);

tpoints = round(T/K);
xpoints = round(2*L/H);
xgrid = linspace(-L_fd,L_fd, xpoints);
tgrid = linspace(0,T,tpoints);




lambdagrid = sigma0*exp(lam*tgrid);

%% Produce sigma values.

figure;
hold on;

% Plot the solution at sigma0
u0 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u0 = u0+sigma0^j*y(1,:);
end
plot(xgrid,u0,'-k','LineWidth',2);

% Plot the solution at sigma1
u1 = zeros(1,length(xgrid));
for j = 0:num_homological_equations
    y = P{j+1}.fun(xgrid);
    u1 = u1+sigma1^j*y(1,:);
end
plot(xgrid,u1,'-g','LineWidth',2);

%% Time evolution

% Set up grid
T = (log(sigma1)-log(sigma0))/3;
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);

U_n = u0;
U_o = u0;
tol = 1e-8;
p = struct;

temp_L = zeros(2,1);
temp_R = zeros(2,1);
for j = 1:length(tgrid)

    
    % Update the boundary values.
    sigman = lambdagrid(j);
    
    bc_L = 0;
    bc_R = 0;
    % Produce the boundary conditions.
    for i = 0:num_homological_equations
        temp_L = P{i+1}.fun(xgrid(1));
        temp_R = P{i+1}.fun(xgrid(end));
        bc_L = bc_L + sigman^i*temp_L(1,:);
        bc_R = bc_R + sigman^i*temp_R(1,:);
%         d_bc_L = d_bc_L + sigman^i*temp_L(2,:);
%         d_bc_R = d_bc_R + sigman^i*temp_R(2,:);
    end
   
    
    % Create boundary functions
    bc_L_fun = @(U_n,U_o,K,H,p)(U_n(1,1)-bc_L(1));
    bc_R_fun = @(U_n,U_o,K,H,p)(U_n(1,end)-bc_R(1));
    bc_L_jac_fun = @(U_n,U_o,K,H,p) 1;
    bc_R_jac_fun = @(U_n,U_o,K,H,p) 1;
    
    % Update finite difference
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    % Plot current state
    j
    clf;
    hold on;
    plot(xgrid,u1,'-g','LineWidth',2);
    plot(xgrid,u0,'-k','LineWidth',2);
    plot(xgrid,U_n,'--r','LineWidth',2);
    drawnow;
    pause(0.1);
    U_o = U_n;

end

diff = max(abs(u1-U_n))


% Theta= 1.8 converges.






