function cf_alpha = convert_from_mu_to_alpha(cf,rho_alpha,alpha_left,alpha_right,mu_left, ...
    mu_right,L_rho_alpha,D_rho_alpha,eta_rho_alpha)

% rigorous enclosure of $\pi$.
pie = iv('pi');

% Theta to parametrize the stadium
theta = linspace(0,sup(2*pie),1001);
theta = iv(theta(1:end-1),theta(2:end));

% The stadium in the variable $\tilde \alpha$.
alpha_tilde = (rho_alpha/2)*exp(1i*theta)+1./((2*rho_alpha)*exp(1i*theta));

% $\alpha$ corresponding to $\tilde \alpha$ on the stadium
alpha = (alpha_left+alpha_right)/2+((alpha_left-alpha_right)/2)*alpha_tilde;

% the variable mu corresponding to $\tilde \alpha$ on the stadium
mu = (-2+2*sqrt(2+alpha.^2))./(1+alpha.^2);

% the variable $\tilde mu$ corresponding to $\tilde \alpha$ on the stadium
mu_tilde = (2*mu-(mu_left+mu_right))/(mu_right-mu_left);

% initialize the Chebyshev polynomial interpolant (in variable mu_tilde)
% evaluated on the stadium (the stadium corresponding to alpha_tilde).
p_eval = cf(1)*iv(ones(size(mu_tilde)))+ cf(2)*mu_tilde;

% Evaluate the rest of the Chebyshev polynomial interpolant (in variable mu_tilde)
% evaluated on the stadium (the stadium corresponding to alpha_tilde).
T_old_old = iv(ones(size(mu_tilde))); % Chebyshev polynomial from two indices back
T_old = mu_tilde; % Chebyshev polynomial from one index back
for j = 3:length(cf)
    T_cheb = 2*mu_tilde.*T_old+T_old_old; % Chebyshev polynomial for the current index
    p_eval = p_eval + cf(j)*T_cheb; % update the interpolant
    T_old_old = T_old; % update Chebyshev polynomial from two indices ago
    T_old = T_cheb; % update Chebyshev polynomial from one index ago
end

M_rho_alpha = iv(maxi(sup(abs(p_eval)))); % Bound on the modulus of the function on the stadium

% Find out how many interpolation nodes we need in the varaible $\tilde
% \alpha$ in order to have small error.
error_alpha = 1;
N_alpha = 3;
tol = iv(1e-16);
while error_alpha > inf(tol)
    N_alpha = N_alpha + 1;
    error_alpha = M_rho_alpha*L_rho_alpha/(pie*D_rho_alpha*sinh(eta_rho_alpha*N_alpha+1));
end

%{
    Now evaluate the function on the interpolation nodes.
%}

% Determine the interpolation nodes
theta = (2*iv(0:1:N_alpha-1)+1)*pie/(2*N_alpha);
alpha_tilde = cos(theta); % nodes in $\tilde \alpha$
alpha = (alpha_left+alpha_right)/2+((alpha_left-alpha_right)/2)*alpha_tilde; % nodes in $\alpha$
mu = (-2+2*sqrt(2+alpha.^2))./(1+alpha.^2); %nodes in $\mu$

mu_tilde = (2*mu-(mu_left+mu_right))/(mu_right-mu_left); % nodes in $\tilde \mu$

% Matrix for getting Chebyshev interpolant coefficents
Id2 = (iv(2)/N_alpha)*speye(N_alpha);
Id2(1,1) = Id2(1,1)/iv(2);
Tcf = Id2*cos(theta.'*(0:1:N_alpha-1)).';
    
% Matrix for evaluating polynomial in $\tilde \mu$ 
theta = acos(mu_tilde);
T = cos(theta.'*(0:1:length(cf)-1));

% Evaluation of polynomial in $\tilde \mu$ at the nodes for $\tilde \alpha$
fun_eval = T*cf;

cf_alpha = Tcf*fun_eval; % coefficients for the Chebyshev polynomial in alpha

cf_alpha(1) = cf_alpha(1)+iv(-tol,tol); % interpolation error

cf_alpha = clip_tail(cf_alpha); % clip the tail












