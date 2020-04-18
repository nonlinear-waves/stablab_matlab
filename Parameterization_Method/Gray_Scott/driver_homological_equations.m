
%% Nagumo - reaction diffusion equation ( U_t = U_{xx} - U + U^3 )

%{
Numericaly solves for the unstable manifold solution to Nagumo's equation
 associated with the unstable traveling wave solution, $\sqrt{2}\sech(x)$. 
%}

addpath(strcat(pwd,'/../bin/'));
addpath(strcat(pwd,'/../bin/dmsuite'));
clc; clear all; close all; beep off; curr_dir = cd;

%% User defined parameters

% System variables
p.gamma = 1/9; % 0 < gamma < 2/9 

% Spatial infinity (chosen so the $|p_0(L)| \leq 1e-8$, where $p_0$ is the profile.
L = 8;

% Scale of the eigenfunction. Chosen so that the last
% homological equation we solve for is on the order of TOL or smaller.
scl = 1.075; 

% The tolerance we request of the bvp solver.
TOL = 1e-12;

%% Set up variables

%dependent variables (alpha*gamma = 1)
p.alpha = 1/p.gamma;
p.Q = sqrt(1-9*p.gamma/2);

%% solve for the eigenfunction
profile_fun = @(x)[1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    -3*p.Q*sqrt(p.gamma)*sinh(x/sqrt(p.gamma))./(1+p.Q*cosh(x/sqrt(p.gamma))).^2;...
    3./(1+p.Q*cosh(x/sqrt(p.gamma))); ...
    3*p.Q*(1/sqrt(p.gamma))*sinh(x/sqrt(p.gamma))./(1+p.Q*cosh(x/sqrt(p.gamma))).^2];
    
s = struct;
% ODE the eigenfunction satisfies
eigfun_ode = @(x,y,lambda)([A_explicit(x,lambda,s,p)*y(1:4); ...
    A_explicit(x+L,lambda,s,p)*y(5:8)]);

% boundary conditions the eigenfunction satisfies
fun_bc = @(ya,yb,lambda)(bc_eigfun(ya,yb,lambda,p));

% guess for the eigenfunction
guess_fun = @(x)[sech(x);-sech(x)*tanh(x);sech(x);-sech(x)*tanh(x)];
guess = @(x)[guess_fun(x);guess_fun(x+L)];

% guess for the eigenvalue
lambda0 = 1;

% solve for the eigenfunction and eigenvalue
solinit = bvpinit(linspace(-L,0,30),guess,lambda0);
options = bvpset('RelTol',TOL,'AbsTol',TOL,'Nmax', 20000);
eigfun_sol = bvp5c(eigfun_ode,fun_bc,solinit,options);
lam = eigfun_sol.parameters; % eigenvalue of linearized PDE

%% Find the homological equations, P{n}.fun

% Unstable wave solution 
fprintf('n=0, Loading function\n\n');
P{1}.fun = profile_fun;

% Eigenfunction
fprintf('n=1, Loading function\n\n');
P{2}.fun = @(x) scl*eig_fun(x,eigfun_sol).';

%
% Solve for the solutions to the homological equations
%


% degree of Chebyshev interpolant of the solution of the homological
% equation
N = 1001;

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




options = bvpset('RelTol', TOL, 'AbsTol', TOL,'Nmax', 20000);
cnt = 0;
n = 1;
p_inf= TOL+1;
while p_inf > TOL
    
    n = n+1;
    fprintf(strcat('n=', num2str(n), ', Generating function\n\n'));
    
    % ODE for the homological equations
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam,p));
    
    % projective boundary conditions for the homological equations
    bc_handle = @(ya,yb)(bc(ya,yb,n,lam,p)); 
    solinit = bvpinit([-L,L],P{n}.fun);
    sol = bvp5c(ode_handle,bc_handle,solinit);
            
    % record data
    P{n+1}.mu_L = sqrt(3*n+1);
    P{n+1}.mu_R = -sqrt(3*n+1);
    P{n+1}.sol = sol;
    
    % Get function handle for Chebyshev interpolant
    cf = Tcf*deval(sol,x).';
    P{n+1}.fun = @(x)(interp_cheby_fun(cf,x,N,a,b));
    
    % Get the max norm of the solution to the homological equation just
    % solved.
    p_inf = max(norm(deval(sol,x)))

end

%% Plot the homological equations
fprintf('Preparing graphs\n\n');
x = linspace(-L,L,1000);

hold on;
for j = 1:length(P)
    y = P{j}.fun(x);
    plot(x,y,'LineWidth',2)
    nrm = max(max(abs(y)));
end

% Record data and save it.
data.P = P;
data.L = L;
data.p = p;
data.lam = lam;
data.scl = scl;
data.TOL = TOL;
cd('../data');
save('P_Gray_Scott');
cd(curr_dir);


