
%{
    This is the main file. 
%}

clc; close all; clear all; beep off; curr_dir = cd; t_start = tic;

% -------------------------------------------------------------------------
% Controls
% -------------------------------------------------------------------------

% left bound on the interval in the parameter nu
nu_left = '0.28';

% right bound on the interval in the parameter nu
nu_right = '0.281';

% The radius in which we seek to verify a unique solution to the BVP
% problems about our initial guess z_bar
r_newton = iv(1e-6);

% The scalar we use for the eigenvalues in the parameterization of the
% stable manifold at positive infinity.
scale_R = 0.5;

% The approximate choice for $\Delta x$ to solve the ODE forward in the
% middle manifolds
approx_delx = 0.05;

plot_it = 'on'; % plot the profile or not

% -------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------

% file name to save
file_name = strrep(['proof_',nu_left,'_',nu_right],'.','p');

nu_left = iv(nu_left); % make left nu value an interval
nu_right = iv(nu_right); % make right nu value an interval

min_x_R = 4; % Require the middle interval to be at least 4 unites wide

% The maximum degree of the Chebyshev polynomials in the parameter mu
M = 20;

% -------------------------------------------------------------------------
% Initialize arrays for the Newton-Kantorovich argument
% -------------------------------------------------------------------------

% determine how wide the middle interval needs to be using nonrigorous
% numerics
[x_R,theta_0,theta_R,r_R] = get_x_R(sup(nu_right),scale_R,min_x_R);
% x_R - the right bound on the middle manifold part
% theta_0 - the parameter for the left manifold
% theta_R - the angle parameter for the right manifold (theta_R is called
% tau in the paper)
% r_R - the radius parameter for the (2D) right manifold

% the number of nodes in the middle manifold part
num_nodes = ceil(x_R/approx_delx)+1;

% The number of equations in the Newton-Kantorovich argument. The ODE has 4
% components, 2 for the profile and 2 for the Ricatti equation for a total
% of 4. This gives 4*num_nodes variables, plus theta_L and theta_R for
% 4n+2 total variables. That is, we will be using the Newton-Kantorovich
% argument to find the zeros of a function $F:\mathbb{R}^{2n+2}\to
% \mathbb{R}^{2n+2}$, whose zeros correspond to the existence of a global
% solution to the BVP problem.
Nequ = 4*num_nodes+2;

% initialize z_bar and F and its zeros, in the parameter mu
z_bar_mu = zeros(Nequ,M); % The guess at the zeros of F
r_R_mu = zeros(M,1); % the radius of the parameter at + infty
F_mu = iv(zeros(Nequ,M)); % The function F, whose zeros we seek
DF_mu = iv(zeros(Nequ,Nequ,M)); % The first derivatives of F
pre_z_bar_mu = zeros(size(z_bar_mu)); % The vector that will be used to get z_bar_mu
pre_r_R_mu = zeros(size(r_R_mu)); % The vector that will be used to get r_R_mu
% fprintf("\nError from right manifold: %4.4g\n ",R_mani_error);

% -------------------------------------------------------------------------
% Make the guess, z_bar, of the zeros of F
% -------------------------------------------------------------------------

theta = (2*(0:1:M-1)+1)*pi/(2*M); % The values of theta for interpolation
mu_tilde = cos(theta); % the values of \tilde mu for interpolation
mu_right = (-1/(2*nu_left))*(1-sqrt(1+4*nu_left));  % right bound on mu
mu_left = (-1/(2*nu_right))*(1-sqrt(1+4*nu_right)); % left bound on mu
mu = mid((mu_left+mu_right)/2+(mu_right-mu_left)*mu_tilde/2); % values of mu for interpolation
nu = (1-mu)./mu.^2; % values of nu for interpolation

% Use double arithmetic to determine the needed values of mu
for j = 1:M

    % determine theta_R, r_R, and the ODE solution for the BVP problem at
    % the interpolatio nodes. 
    [theta_R,r_R,sol_u] = fixed_x_R(nu(j),scale_R, ...
    x_R,theta_0,theta_R,r_R);

    % evlaute the BVP solution at the spatial nodes
    mid_domain = linspace(0,x_R,num_nodes);
    temp = deval(sol_u,mid_domain);
   
    % record the guess, $\bar z$, of the zeros of $F$
    pre_z_bar_mu(:,j) = [theta_0;reshape(temp,size(temp,1)*size(temp,2),1);theta_R];

    % record the radius of the parameter for evaluation of the right
    % manifold
    pre_r_R_mu(j) = r_R;

    % Plot the profile and the solution of the Ricatti equation, if desired
    if strcmp(plot_it,'on')
        figure
        hold on;s
        plot(sol_u.x,sol_u.y,'-k','LineWidth',2);
        plot(mid_domain,temp,'.r','MarkerSize',18);
        h = xlabel('x');
        set(h,'FontSize',18);
        h = ylabel('\phi, \phi_x, u, u_x');
        set(h,'FontSize',18);
        h = gca;
        set(h,'FontSize',18);
    end

    return

end

% Now obtain z_bar with entries that are Chebyshev polynomials in the
% variable mu
Tcf = mid(get_Tcf(M));
z_bar_mu(1,1) = pre_z_bar_mu(1); % Keep theta_0 constant
for j = 2:size(pre_z_bar_mu,1)
    z_bar_mu(j,1:M) = Tcf*squeeze(pre_z_bar_mu(j,1:M)).';
end
z_bar_mu = iv(z_bar_mu);

% Now obtain r_R_mu as a Chebyshev polynomial in mu
r_R_mu = Tcf*pre_r_R_mu(1:M);
r_R_mu = iv(r_R_mu);

% -------------------------------------------------------------------------
% Rigorously compute the left manifold
% -------------------------------------------------------------------------

dl = left_manifold(nu_left,nu_right,z_bar_mu(1,1));

dl_wide = left_manifold(nu_left,nu_right,z_bar_mu(1,1)+iv(-r_newton,r_newton));

% -------------------------------------------------------------------------
% Fill in F and its derivatives for the left intervals
% -------------------------------------------------------------------------

% first equation of F: Phi(theta)_1 - phi_0 = Phi(theta)_1-z_bar(2)
temp = poly_add(dl.cf_phi,-z_bar_mu(2,:));
F_mu(1,1:length(temp)) = temp; 

temp = dl.cf_phi_theta;
DF_mu(1,1,1:length(temp)) = temp; % The derivative of the first equation with resepct to theta_0
DF_mu(1,2,1) = -1; % The derivative of the first equation with respect to phi_0

% second equation of F: Phi(theta)_2-w_0 = Phi(theta)_2 - z_bar(3)
temp = poly_add(dl.cf_w,-z_bar_mu(3,:));
F_mu(2,1:length(temp)) =  temp;

temp = dl.cf_w_theta;
DF_mu(2,1,1:length(temp)) = temp; % The derivative of the second equation with resepct to theta_0
DF_mu(2,3,1) = -1; % The derivative of the first equation with respect to w_0

% bound on the left manifold error
L_mani_error = max(sup(abs(F_mu(1:2,1))));

% -------------------------------------------------------------------------
% Rigorously compute the right manifold and fill in F and its derivatives
% -------------------------------------------------------------------------

% compute the right manifold
disp('Computing the right manifold.');
dr = right_manifold(nu_left,nu_right,scale_R);

% transform the parameter nu to the parameter alpha
alpha_right = sqrt(iv(sup(4*nu_right-1)));
if inf(4*nu_left-1)<0
    alpha_left = iv(0);
else
    alpha_left = sqrt(4*nu_left-1);
end

% prepare to do analytic interpolation to convert from mu to alpha
rho_alpha = iv('0.9')*(9/(5*(alpha_right-alpha_left))+ ...
    sqrt(81/(25*(alpha_right-alpha_left)^2)+1));
L_rho_alpha = iv('pi')*sqrt(rho_alpha^2+1/rho_alpha^2);
D_rho_alpha = (rho_alpha+1/rho_alpha)/2-1;
eta_rho_alpha = log(rho_alpha);

% convert theta_R to a polynomial in alpha
theta_R_alpha = convert_from_mu_to_alpha(squeeze(z_bar_mu(end,:)).',rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);

% convert r_R to a polynomial in alpha
r_R_alpha = convert_from_mu_to_alpha(r_R_mu,rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);

% evaluate the right manifold
disp('Evaluating the right manifold.');

rm_eval = eval_rm(theta_R_alpha,r_R_alpha,dr);

disp('Done evaluating the right manifold.');


% RHS phi
temp = convert_from_mu_to_alpha(squeeze(z_bar_mu(end-4,:)).',rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);
temp = poly_add(rm_eval.y_phi,-temp);
F_mu(end-3,1:length(temp)) = temp;

% RHS w
temp = convert_from_mu_to_alpha(squeeze(z_bar_mu(end-3,:)).',rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);
temp = poly_add(rm_eval.y_w,-temp);
F_mu(end-2,1:length(temp)) = temp;

% RHS u
temp = convert_from_mu_to_alpha(squeeze(z_bar_mu(end-2,:)).',rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);
temp = poly_add(rm_eval.y_u,-temp);
F_mu(end-1,1:length(temp)) = temp;

% RHS z
temp = convert_from_mu_to_alpha(squeeze(z_bar_mu(end-1,:)).',rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha);
temp = poly_add(rm_eval.y_z,-temp);
F_mu(end,1:length(temp)) = temp;

% error from right manifold
R_mani_error = max(sup(abs(F_mu(end-3:end,1))));

error_check = max(max(sup(abs(F_mu))));
fprintf('\nRight manifold error: %4.4g\n\n',error_check);

% evaluate the second derivative of the right manifold on a ball or radius
% r_newton
theta_R_alpha_wide = theta_R_alpha;
theta_R_alpha_wide(1) = theta_R_alpha_wide(1)+iv(-r_newton,r_newton); 
rm_wide_eval = eval_rm(theta_R_alpha_wide,r_R_alpha,dr);

% update the Jacobian
for j = 1:4
    DF_mu(end-j+1,end-j,1) = -1;
end

% update the Jacobian
DF_mu(end-3,end,1:length(rm_eval.y_phi_theta)) = rm_eval.y_phi_theta;
DF_mu(end-2,end,1:length(rm_eval.y_w_theta)) = rm_eval.y_w_theta;
DF_mu(end-1,end,1:length(rm_eval.y_u_theta)) = rm_eval.y_u_theta;
DF_mu(end,end,1:length(rm_eval.y_z_theta)) = rm_eval.y_z_theta;

% -------------------------------------------------------------------------
% Input the middle manifolds into F and its derivatives
% -------------------------------------------------------------------------

num_mid_equations = length(mid_domain)-1;
I4 = iv(eye(4)); % 4x4 identity
B = 0; % bound on the second derivative
r_int = iv(-r_newton,r_newton); % interval width on which we bound the second derivative
states = zeros(num_mid_equations,3);
dmm(num_mid_equations).dm = struct();

disp('Computing the middle manifold.');
for j = 1:num_mid_equations

    j

    delx = iv(mid_domain(j+1))-mid_domain(j);
    
    % get the initial conditions
    phi0 = squeeze(z_bar_mu(4*j-2,:));
    w0 = squeeze(z_bar_mu(4*j-1,:));
    u0 = squeeze(z_bar_mu(4*j,:));
    z0 = squeeze(z_bar_mu(4*j+1,:));

    % compute the middle manifold on a step of size Delta x
    dm = mid_section_solver(phi0,w0,u0,z0,mu_left,mu_right,delx,r_newton);
    
    % update the states vector to determine if there are any bound states
    states(j,:)=[dm.u0_vs_0 dm.u_dx_vs_0 dm.z_dx_vs_0];

    dmm(j).dm = dm;
    
    ab = squeeze(dmm(j).dm.poly_D2Psi_n(:,:,:,:));
    if sum(sum(sum(sum(isnan(ab))))) > 0
       return 
    end
    
    % phi
    temp = poly_add(dm.poly_phi,-z_bar_mu(4*j+2,:));
    F_mu(4*j-1,1:length(temp)) = temp;
    
    % w
    temp = poly_add(dm.poly_w,-z_bar_mu(4*j+3,:));
    F_mu(4*j,1:length(temp)) = temp;
    
    % u
    temp = poly_add(dm.poly_u,-z_bar_mu(4*j+4,:));
    F_mu(4*j+1,1:length(temp)) = temp;
    
    % z
    temp = poly_add(dm.poly_z,-z_bar_mu(4*j+5,:));
    F_mu(4*j+2,1:length(temp)) = temp;
    
    mid_mani_error = max(sup(abs(F_mu(4*j-1:4*j+2,1))));
    fprintf("\n\nError from the middle manifold: %4.4g\n ",mid_mani_error);

    % update the Jacobian
    DF_mu(4*j-1:4*j+2,4*j-2:4*j+1,1:size(dm.poly_DPsi,3)) = dm.poly_DPsi;
    
    DF_mu(4*j-1:4*j+2,4*j+2:4*j+5,1) = -I4;
    
end

% -------------------------------------------------------------------------
% Determine if there is more than one conjugate point
% -------------------------------------------------------------------------

% Set state_now to be the sign of u at the first node.
% If the sign is undetermined, pop out an error. Note that, in below, state_now is
% always either 1 or -1, and never set to be zero.
state_now=states(1,1);
if state_now==0
    error('cannot validate crossing')
end
% u_crossing counts the number of conjugate points.
% Initialize it to be zero.
u_crossings = 0;
for j = 2:num_mid_equations
    % Check if state at the jth node is the same as the current state
    % state_now.
    if states(j,1)==state_now
        % If so, validate u is sign-definite on the interval between
        % the (j-1)th node and the jth node by either validating that u is
        % monotone on the interval or u itself is sign-definite. If the
        % validation passes, continue to the next interval, otherwise pop
        % out an error.
        if states(j-1,3)~=0||states(j-1,2)==state_now
            continue
        else
            error('cannot validate crossing')
        end
        % If the sign of u is undetermined at the jth node, we must
        % validate u is monotone on the the interval between the (j-1)th node
        % and the jth node, otherwise we pop out an error.
    elseif states(j,1)==0
        if states(j-1,3)==0
            error('cannot validate crossing')
        end
        % If u is sign-definite at the jth node and the sign of u is different
        % from current state state_now, validate u is monotone on the the
        % interval between the (j-1)th node and the jth node and then add 1 to
        % u_crossing and update the current state to state(j,1), otherwise pop out an error.
    else
        if states(j-1,3)~=0
            u_crossings=u_crossings+1;
            state_now=states(j,1);
        else
            error('cannot validate crossing')
        end
    end
end

if ~(u_crossings == 1)
    error('There was not exactly one u_crossing in the middle manifold.');
end

% -------------------------------------------------------------------------
% Newton-Kantorovich argument
% -------------------------------------------------------------------------

% the number of x-intervals to use in the computation
intervals = 50;

% -------------------------------------------------------------------------
% form mu_theta and alpha_theta
% -------------------------------------------------------------------------

mu_theta = linspace(0,sup(iv('pi')),intervals+1);
theta_boundaries = iv(mu_theta);
mu_theta = iv(mu_theta(1:end-1),mu_theta(2:end));

mu_boundaries = (mu_left+mu_right)/2+(mu_right-mu_left)*cos(theta_boundaries)/2;

alpha_boundaries = sqrt(4-4*mu_boundaries-mu_boundaries.^2)./mu_boundaries;

alpha_tilde_boundaries = (2*alpha_boundaries-(alpha_left+alpha_right))/(alpha_right-alpha_left);
alpha_tilde_boundaries(1) = iv(-1,sup(alpha_tilde_boundaries(1)));
alpha_tilde_boundaries(end) = iv(inf(alpha_tilde_boundaries(end)),1);

theta_alpha_boundaries = acos(alpha_tilde_boundaries);
alpha_theta = iv(theta_alpha_boundaries(1:end-1),theta_alpha_boundaries(2:end));

% -------------------------------------------------------------------------
% Convert from Chebyshev coefficients to x-intervals
% -------------------------------------------------------------------------

%
% Set up the transformation
%

% Get the dimensions of the Jacobian DF_mu
dim = size(DF_mu,1); % Dimension of the Jacobian
N_DF_cheb = size(DF_mu,3); % number of Chebyshev coefficients

% Form T, the transformation matrix between Chebyshev coefficients and
% x-intervals
half = iv(1)/2; % 1/2
pie = iv('pi'); % pi

T_DF_mu = cos(mu_theta.'*(0:1:N_DF_cheb-1)); % cos(theta*n)
T_DF_alpha = cos(alpha_theta.'*(0:1:N_DF_cheb-1)); % cos(theta*n)

%
% Get DF, the evaluation of DF_mu at the x-intervals
%

DF = iv(zeros(dim,dim,intervals));
% D2F = iv(zeros(dim,dim,dim,intervals));
for j = 1:dim-4
        DF(j,:,:) = (T_DF_mu*squeeze(DF_mu(j,:,:)).').';
end

for j = dim-3:dim
    DF(j,:,:) = (T_DF_alpha*squeeze(DF_mu(j,:,:)).').';
end
       
B = 0;
for j = 1:length(dmm)
    for q = 1:4
        for k = 1:4
            deg = size(dmm(j).dm.poly_D2Psi_n,4);
            T_D2F = cos((0:1:deg-1).'*mu_theta); % cos(theta*n)
            temp = squeeze(dmm(j).dm.poly_D2Psi_n(q,k,:,:))*T_D2F;
            B = maxi([B,maxi(maxi(sup(abs(temp))))]);
        end
    end
end

temp = cos(mu_theta.'*(0:1:length(dl_wide.cf_phi_theta_theta)-1))*dl_wide.cf_phi_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

temp = cos(mu_theta.'*(0:1:length(dl_wide.cf_w_theta_theta)-1))*dl_wide.cf_w_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

temp = cos(mu_theta.'*(0:1:length(rm_wide_eval.y_phi_theta_theta)-1))*rm_wide_eval.y_phi_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

temp = cos(mu_theta.'*(0:1:length(rm_wide_eval.y_w_theta_theta)-1))*rm_wide_eval.y_w_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

temp = cos(mu_theta.'*(0:1:length(rm_wide_eval.y_u_theta_theta)-1))*rm_wide_eval.y_u_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

temp = cos(mu_theta.'*(0:1:length(rm_wide_eval.y_z_theta_theta)-1))*rm_wide_eval.y_z_theta_theta;
B = maxi([B,maxi(sup(abs(temp)))]);

B = iv(B);

%
% Get F, the evaluation of F_mu at the x-intervals
%

% Form T, the transformation matrix between Chebyshev coefficients and
% x-intervals

N_F_cheb = size(F_mu,2);
T_F = cos(mu_theta.'*(0:1:N_F_cheb-1)); % cos(theta*n)

F = (T_F*squeeze(F_mu(:,:)).').';

% -------------------------------------------------------------------------
% Get bounds for the Newton-Kantorovich argument
% -------------------------------------------------------------------------

J = mid(DF); % midpoint of DF

% Choose the approximate inverse of A_dagger
A = zeros(size(J));
for j = 1:intervals
    A(:,:,j) = squeeze(J(:,:,j))\eye(size(J,1));
end
A = iv(A);

% Choose A_dagger to be the midpoint of DF
A_dagger = iv(J);

% Obtain an upper bound Z0 on Id-A*A_dagger on the x-intervals
Id = iv(eye(dim));
max_Z0 = 0;
for j = 1:intervals
   max_Z0 = maxi(max_Z0,sup(norm(Id-squeeze(A(:,:,j)) ...
       *squeeze(A_dagger(:,:,j)),2)));
end
Z0 = iv(max_Z0)


% Obtain an upper bound Y0 on A*F on the x-intervals
max_Y0 = 0;
for j = 1:intervals
   max_Y0 = maxi(max_Y0,sup(norm(squeeze(A(:,:,j))*F(:,j),2)));
end
Y0 = iv(max_Y0)

% Obtain an upper bound Z1 on A*(A_dagger-DF) on the x-intervals
max_Z1 = 0;
for j = 1:intervals
   max_Z1 = maxi(max_Z1,sup(norm(squeeze(A(:,:,j))* ...
       (squeeze(A_dagger(:,:,j))-squeeze(DF(:,:,j))),2)));
end
Z1 = iv(max_Z1)

max_Z2r = 0;
for j = 1:intervals
    max_Z2r = maxi(max_Z2r,sup(norm(squeeze(A(:,:,j)),'fro')));
end
Z2r = (2*B)*dim*iv(max_Z2r) % There are up to 4 variables that have non-zero derivatives

temp = (1-Z0-Z1)^2-4*Z2r*Y0;
if inf(temp)< 0
   error('fatal failure'); 
end
r_NK = iv(min(inf(iv(0.9*(mid(((1-Z0-Z1)+sqrt(inf(temp)))/(2*Z2r))))),inf(r_newton)));

if sup(r_NK) > inf(r_newton)
       error('fatal failure'); 
end
 
newton_kantorovich_poly = Z2r*r_NK^2-(1-Z0-Z1)*r_NK+Y0;

if sup(newton_kantorovich_poly) >= 0
   success = 0;
else
   success = 1;
end

if success == 1
   disp('Congratulations! The proof succeeded.'); 
end

total_run_time = toc(t_start);

% -------------------------------------------------------------------------
% Save the data
% -------------------------------------------------------------------------

cd('../data');
save(file_name);
cd('../code');


