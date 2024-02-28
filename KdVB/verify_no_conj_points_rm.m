function verify_no_conj_points_rm(r,dr,x_R,alpha_left,alpha_right,theta_points)


% clc;

% cd('../data');
% ld = load('nu_left0_52nu_right0_53.mat');
% cd('../code');


% theta_points = 101;
% r = ld.r_R_alpha;
% dr = ld.dr;
% x_R = ld.x_R;
% alpha_left = ld.alpha_left;
% alpha_right = ld.alpha_right;


%
% Rigorously represent $1/(1+alpha^2)$ in terms of $\tilde alpha$.
%

% alpha = (alpha_L+alpha_R)/2+(alpha_L-alpha_R)*tilde_alpha/2
alpha = [(alpha_left+alpha_right)/2,(alpha_left-alpha_right)/2];

% the square of alpha as a polynomial in the variable tilde_alpha
alpha_squared = poly_mult(alpha,alpha);

% 1+alpha^2 as a polynomial in the variable tilde_alpha
denom = poly_add(iv([1;0]),alpha_squared); % 1 + alpha^2

U0 = denom(1);
del_U = iv(zeros(size(denom)));
del_U(2:end) = denom(2:end);

% Taylor expansion of 1/(U0+del_U) at U0 in the variable del_U
term_pow = 8;
del_U_pow_j = iv([1;0]);
term = 1/U0;
for j = 1:term_pow

    % del_U^j
    del_U_pow_j = poly_mult(del_U,del_U_pow_j);

    % derivative of 1/(U0+del_U) with respect to del_U evaluated at U0 
    temp = ((-1)^j/U0^(j+1))*del_U_pow_j;

    % partial sum of the Taylor expansion of 1/(U0+del_U)
    term = poly_add(term,temp);
end

% Taylor remainder term
remainder = (-1)^(term_pow+1)*poly_mult(del_U_pow_j,del_U);

% Rigorous enclosure of 1/(U0+del_U) 
term = poly_add(term,remainder);

% get the real and imaginary part of $\mu_{+}$
real_mu_plus = clip_tail(-2*term);

%{

Let theta_+^0 be the parameter for the right manifold in the application of
the Newton-Kantorovich Theorem. Then theta_+^0 = r_0*(cos(theta)+i*sin(theta).
Now theta_+(x) = theta_+^0*exp(mu_+*(x-x_R)).  Thus

theta_+(x) = r(x)*(cos(theta(x))+1i*sin(theta(x)))

where

r(x) = |theta_+^0|*exp(real(mu_+)*(x-x_R))

and

theta(x) = theta_+^0 + imag(mu_+)(x-x_R)

%}

% let theta range over [0,2*pi]
theta_inf_vals = linspace(0,sup(2*iv('pi')),theta_points);
theta_inf_vals = iv(theta_inf_vals(1:end-1),theta_inf_vals(2:end));

for j = 1:length(theta_inf_vals)

    j

    theta_inf = theta_inf_vals(j);
        
    % find how big the radius can be for $\theta_{\pm}$. 
    
    % real_mu_plus is a polynomial in the variable tilde_alpha. Bound it above
    bound_on_real_mu_plus = real_mu_plus(1) + sum(abs(real_mu_plus(2:end)));
    if sup(bound_on_real_mu_plus) >= 0
        error('failed to verify');
    end
    % get bound on e^{mu*x} for $x\in [x_R,\infty)$
    temp = (bound_on_real_mu_plus)*(x_R);
    exp_real_mu_plus_x_sup = sup(exp(temp));
    
    % enclose the radius
    r_inf = poly_mult(r,iv(0,exp_real_mu_plus_x_sup));
    
    % evaluate the right manifold
    d_temp = eval_rm(theta_inf,r_inf,dr);
    
    % bound $u(x)$ from below for $x\in[x_R,+\infty)$
    condition = inf(d_temp.y_u(1)-sum(abs(d_temp.y_u(2:end)))) ;
    
    if isnan(condition)
        error('Failed to verify that there are no conjugate points in the right manifold');
    end
    
    if condition < 0
        error('Failed to verify that there are no conjugate points in the right manifold');
    end
        
end












