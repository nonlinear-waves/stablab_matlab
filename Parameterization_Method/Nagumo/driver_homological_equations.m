
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

% The tolerance we request of the bvp solver.
TOL = 1e-12;

%% Set up variables

lam = 3; % eigenvalue of linearized PDE

%% Find the homological equations, P{n}.fun

% Unstable wave solution 
fprintf('n=0, Loading function\n\n');
P{1}.fun = @(x)([sqrt(2)*sech(x); ...
                sqrt(2)*(-tanh(x).*sech(x))]);

% Eigenfunction
fprintf('n=1, Loading function\n\n');
P{2}.fun = @(x) scl*eig_fun(x);

% % % %
% % % % Solve for the solutions to the homological equations
% % % %
% % % 
% % % options = bvpset('RelTol', 1e-12, 'AbsTol', 1e-12,'Nmax', 20000);
% % % cnt = 0;
% % % for n = n_vals
% % %     
% % %     fprintf(strcat('n=', num2str(n), ', Generating function\n\n'));
% % %     
% % %     ode_handle = @(x,y)(ode_fun(x,y,P,n,lam));
% % %     
% % %     bc_handle = @(ya,yb)(bc_fun(ya,yb,n)); 
% % %     solinit = bvpinit([-L,L],P{n}.fun);
% % %     sol = bvp5c(ode_handle,bc_handle,solinit);
% % %             
% % %     P{n+1}.mu_L = sqrt(3*n+1);
% % %     P{n+1}.mu_R = -sqrt(3*n+1);
% % %     P{n+1}.sol = sol;
% % %     
% % %     % degree of polynomial
% % %     N = 240;
% % %     
% % %     a = -L;
% % %     b = L;
% % % 
% % %     % Chebyshev nodes
% % %     theta = ((0:1:N-1)+0.5)*pi/N;
% % %     x0 = cos(theta); % in [-1,1]
% % %     x = 0.5*(a+b)+0.5*(a-b)*x0; % in [a,b]
% % % 
% % %     % Transformation to get Chebyshev coefficients
% % %     Id2 = (2/N)*speye(N);
% % %     Id2(1,1) = Id2(1,1)/2;
% % %     Tcf = Id2*cos(theta.'*(0:1:N-1)).';
% % %     cf = Tcf*deval(sol,x).';
% % %     
% % %     P{n+1}.fun = @(x)(interp_cheby_fun(cf,x,n,N,a,b));
% % %     
% % %     nrm = max(norm(deval(sol,x)))
% % % 
% % % end
% % % 
% % % % %% Plot the homological equations
% % % % fprintf('Preparing graphs\n\n');
% % % % x = linspace(-L,L,1000);
% % % % 
% % % % hold on;
% % % % for j = 1:num_homological_equations
% % % %     y = P{j}.fun(x);
% % % %     plot(x,y,'LineWidth',2)
% % % %     nrm = max(max(abs(y)));
% % % % end
% % % % 
% % % % data.P = P;
% % % % data.L = L;
% % % % data.N = N;

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
    ode_handle = @(x,y)(ode_fun(x,y,P,n,lam));
    
    % projective boundary conditions for the homological equations
    bc_handle = @(ya,yb)(bc_fun(ya,yb,n)); 
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
data.lam = lam;
data.scl = scl;
data.TOL = TOL;
cd('../data');
save('P_Nagumo','data');
cd(curr_dir);



