
%% Nagumo - reaction diffusion equation ( U_t = U_{xx} - U + U^3 )

%{
Numericaly solves for the unstable manifold solution to Nagumo's equation
 associated with the unstable traveling wave solution, $\sqrt{2}\sech(x)$. 
%}

addpath(strcat(pwd,'/../bin/'));
addpath(strcat(pwd,'/../bin/dmsuite'));
clc; clear all; close all; beep off; 
curr_dir = cd;

%% User defined parameters

% Scale of the eigenfunction. Chosen so that the last
% homological equation we solve for is on the order of TOL or smaller.
scl = 0.1; 

%% Find the homological equations, P{n}.fun

% Eigenfunction
cd('../data');
ld = load('eig_fun_Schrodinger_a_4_gamma_2');
data = ld.data;
N = data.N;
p = data.p;
cd(curr_dir);




% return

% x = linspace(data.a,data.b,4001);
% plot(x,eval_eigfun(data.sol,x,data.b))


% numerical infinity
L = data.b;

% Unstable wave solution 
fprintf('n=0, Loading function\n\n');

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


% Get function handle for Chebyshev interpolant
cf = Tcf*profile_fun(x).';
P{1}{1}.fun = @(x)(interp_cheby_fun(cf,x,N,a,b));


% Get function handle for Chebyshev interpolant
cf = Tcf*scl*eval_eigfun(data.sol,x,data.b).';
P{2}{1}.fun = @(x)(interp_cheby_fun(cf,x,N,a,b));
P{1}{2}.fun = @(x)conj(P{2}{1}.fun(x));

dom = linspace(a,b,8309);
err = max(max(abs(P{2}{1}.fun(dom)-scl*eval_eigfun(data.sol,dom,data.b))))

% eigenvalue associated with the eigenfunction
lam1 = data.lam;
lam2 = conj(lam1);

options = bvpset('RelTol', 1e-6, 'AbsTol', 1e-8,'Nmax', 20000,'Vectorize', ...
        'off','stats','on');

for diag_number = 2:20
    for n = 0:diag_number
        m = diag_number - n;

%         if n > m
%            P{n+1}{m+1}.fun = @(x)conj(P{m+1}{n+1}.fun(x));
%            continue
%         end
        
        
        fprintf('\nSolving the Homological equations for (n,m) = (%g,%g)\n',n,m);
    
        % ODE for the homological equations
        ode_handle = @(x,y) ode_fun(x,y,n,m,lam1,lam2,P,p);

        % projective boundary conditions for the homological equations
        bc_handle = @(ya,yb) bc(ya,yb,n,m,lam1,lam2,p); 

        guess_ind_1 = max(2,n)-1;
        guess_ind_2 = max(2,m)-1;
        init_guess = @(x) P{guess_ind_1}{guess_ind_2}.fun(x);

        solinit = bvpinit(linspace(-L,L,30),init_guess);

        sol = bvp5c(ode_handle,bc_handle,solinit);


        % Get function handle for Chebyshev interpolant
        cf = Tcf*deval(sol,x).';
        P{n+1}{m+1}.fun = @(x)(interp_cheby_fun(cf,x,N,a,b));
        
%         P{n+1}{m+1}.fun = makeChebychevFast(@(x)deval(sol,x), L*10, L, 10);

        % Get the max norm of the solution to the homological equation just
        % solved.
        p_inf = max(norm(P{n+1}{m+1}.fun(linspace(-L,L,1001))))
        p_end_norm = norm(P{n+1}{m+1}.fun(a))

    end
end

%% Plot the homological equations
% fprintf('Preparing graphs\n\n');
% x = linspace(-L,L,1000);
% 
% hold on;
% for j = 1:length(P)
%     for k = 1:j
%         y = P{j}{k}.fun(x)+P{k}{j}.fun(x);
%         plot(x,real(y),'LineWidth',2)
%         nrm = max(max(abs(y)));
%     end
% end

% Record data and save it.
data.P = P;
data.L = L;
data.p = p;
data.lam = lam1;
data.scl = scl;
cd('../data');
save('P_Schrodinger','data');
cd(curr_dir);


