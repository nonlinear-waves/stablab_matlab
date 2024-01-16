function d = mid_section_solver(phi0,w0,u0,z0,mu_left,mu_right,delta_x,r_NK)
%{
    The input phi0, w0, u0, z0 are resepectively the Chebyshev coefficints
for the intial condition for phi, w = phi_x, u, and z = u_x. The inputs
mu_left and mu_right are the left and right boundary of the interval on
which the Chebyshev interpolants are given. Note that the Chebyshev
interpolation is in the variable $\mu$.
%}

max_N = 32;
err_tol = 1e-16;


start_time = tic;

% Rigorous enclosure of constants
pie = iv('pi');
half = iv(1)/2;

phi_fun = @(mu)eval_cf_z(mu,phi0,mu_left,mu_right);
w_fun = @(mu)eval_cf_z(mu,w0,mu_left,mu_right);
u_fun = @(mu)eval_cf_z(mu,u0,mu_left,mu_right);
z_fun = @(mu)eval_cf_z(mu,z0,mu_left,mu_right);

%check if u0 is sign-definite for mu in [mu_left,mu_right]
NN=101;
mu_a= linspace(mu_left,mu_right,NN);
mu_a = iv(mu_a(1:end-1),mu_a(2:end));
u_i=add_neighborhood(u_fun(mu_a(1)),r_NK);
if in(0,u_i)
    d.u0_vs_0=0;
elseif inf(u_i)>0
    d.u0_vs_0=1;
else
    d.u0_vs_0=-1;
end
i=2;
while d.u0_vs_0~=0 && i<NN
    u_i=add_neighborhood(u_fun(mu_a(i)),r_NK);
    if in(0,u_i)
        d.u0_vs_0=0;
    elseif (inf(u_i)>0)*d.u0_vs_0==-1
        d.u0_vs_0=0;
    end
    i=i+1;
end

%
% radius of the stadium for Chebyshev interpolation
%

%{
     mu_tilde \in [-1,1]
     mu = (mu_left+mu_right)/2+(mu_right-mu_left)*mu_tilde/2
     Require that $\Re(\mu)\leq 1$.
%}
rho_upper_bound = 2*(1-(mu_left+mu_right)/2+sqrt((mu_left-1)*(mu_right-1)))/(mu_right-mu_left);

% radius of the stadium
rho = iv(sup(0.5*rho_upper_bound));

% Check that rho > 1
if inf(rho) <= 1
    error('rho too small');
end

% Check that $Re(mu) < 1$.
test_term = (mu_left+mu_right)/2+(mu_right-mu_left)*(rho+1/rho)/4;
if sup(test_term) >= 1
    error('rho too big');
end

%
% Chebyshev interpolation error
%

% find d how many Chebyshev nodes are needed for the desired error bound
D_rho = (rho+1/rho)/2-1;
L_rho = pie*sqrt(rho^2+1/rho^2);

eta = log(rho);

% form the stadium
theta = linspace(0,sup(2*pie),1000);
theta = iv(theta(1:end-1),theta(2:end));
mu_stadium = (mu_left+mu_right)/2+(mu_right-mu_left)*(rho*exp(1i*theta)+exp(-1i*theta)/rho)/4; % stadium
nu_stadium = (1-mu_stadium)./mu_stadium.^2;

inverse_nu_stadium = 1./nu_stadium;

%
% evaluate initial data of phi, phi_x, u, u_x, and the Jacobian on the stadium
%
Psi=iv(zeros(4,max_N+1,length(mu_stadium)));

Psi(1,1,:) = phi0(1)+phi0(2)*mu_stadium;
Psi(2,1,:) = w0(1)+w0(2)*mu_stadium;
Psi(3,1,:) = u0(1)+u0(2)*mu_stadium;
Psi(4,1,:) = z0(1)+z0(2)*mu_stadium;

Tn = mu_stadium;
Tn_minus = iv(1);
for j = 3:length(phi0)
    T_cheb = 2*mu_stadium.*Tn-Tn_minus; % next Chebyshev poly
    Tn_minus = Tn; % update old Chebyshev poly
    Tn = T_cheb; % update current Chebyshev poly
    Psi(1,1,:) = transpose(squeeze(Psi(1,1,:))) + phi0(j)*T_cheb; % evaluate the Chebyshev poly
    Psi(2,1,:) = transpose(squeeze(Psi(2,1,:))) + w0(j)*T_cheb;
    Psi(3,1,:) = transpose(squeeze(Psi(3,1,:))) + u0(j)*T_cheb;
    Psi(4,1,:) = transpose(squeeze(Psi(4,1,:))) + z0(j)*T_cheb;
end

Psi_n=iv(zeros(4,max_N+1,length(mu_stadium)));

Psi_n(1,1,:) = add_neighborhood(Psi(1,1,:),r_NK);
Psi_n(2,1,:) = add_neighborhood(Psi(2,1,:),r_NK);
Psi_n(3,1,:) = add_neighborhood(Psi(3,1,:),r_NK);
Psi_n(4,1,:) = add_neighborhood(Psi(4,1,:),r_NK);
% Now compute the coefficients of the series representation of
% Phi on the stadium (variable mu)


Psi(1,2,:) = Psi(2,1,:);
Psi(2,2,:) = inverse_nu_stadium.*transpose(squeeze(-Psi(2,1,:))+half*(squeeze(Psi(1,1,:)).^2-1));
Psi(3,2,:) = Psi(4,1,:);
Psi(4,2,:) = half*Psi(2,1,:).*Psi(3,1,:);

Psi_n(1,2,:) = Psi_n(2,1,:);
Psi_n(2,2,:) = inverse_nu_stadium.*transpose(squeeze(-Psi_n(2,1,:))+half*(squeeze(Psi_n(1,1,:)).^2-1));
Psi_n(3,2,:) = Psi_n(4,1,:);
Psi_n(4,2,:) = half*Psi_n(2,1,:).*Psi_n(3,1,:);

for i = 2:max_N
    
    Psi(1,i+1,:) = Psi(2,i,:)/i;
    Psi(2,i+1,:) = inverse_nu_stadium.*(-transpose(squeeze(Psi(2,i,:)))...
        +half*sum(squeeze(Psi(1,1:i,:)).*squeeze(Psi(1,i:-1:1,:))))/i;
    Psi(3,i+1,:) = Psi(4,i,:)/i;
    Psi(4,i+1,:) = sum(squeeze(Psi(2,1:i,:)).*squeeze(Psi(3,i:-1:1,:)))/(2*i);

    Psi_n(1,i+1,:) = Psi_n(2,i,:)/i;
    Psi_n(2,i+1,:) = inverse_nu_stadium.*(-transpose(squeeze(Psi_n(2,i,:)))...
        +half*sum(squeeze(Psi_n(1,1:i,:)).*squeeze(Psi_n(1,i:-1:1,:))))/i;
    Psi_n(3,i+1,:) = Psi_n(4,i,:)/i;
    Psi_n(4,i+1,:) = sum(squeeze(Psi_n(2,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:)))/(2*i);
    
end


%initializing the jacobian matrix of Psi with respect to Psi0
DPsi = iv(zeros(4,4,max_N+1,length(mu_stadium)));
DPsi(1,1,1,:)=iv(1);
DPsi(2,2,1,:)=iv(1);
DPsi(3,3,1,:)=iv(1);
DPsi(4,4,1,:)=iv(1);

DPsi_n = iv(zeros(4,4,max_N+1,length(mu_stadium)));
DPsi_n(1,1,1,:)=iv(1);
DPsi_n(2,2,1,:)=iv(1);
DPsi_n(3,3,1,:)=iv(1);
DPsi_n(4,4,1,:)=iv(1);

%computing the jacobian matrix on the stadium
for i=1:max_N
    for j=1:4
        
        DPsi(1,j,i+1,:)=DPsi(2,j,i,:)/i;
        DPsi(2,j,i+1,:)=inverse_nu_stadium.*(-transpose(squeeze(DPsi(2,j,i,:)))...
            +sum(squeeze(Psi(1,1:i,:)).*squeeze(DPsi(1,j,i:-1:1,:))))/i;
        DPsi(3,j,i+1,:)=DPsi(4,j,i,:)/i;
        DPsi(4,j,i+1,:)=sum(squeeze(DPsi(2,j,1:i,:)).*squeeze(Psi(3,i:-1:1,:))...
            +squeeze(Psi(2,1:i,:)).*squeeze(DPsi(3,j,i:-1:1,:)))/(2*i);

        DPsi_n(1,j,i+1,:)=DPsi_n(2,j,i,:)/i;
        DPsi_n(2,j,i+1,:)=inverse_nu_stadium.*(-transpose(squeeze(DPsi_n(2,j,i,:)))...
            +sum(squeeze(Psi_n(1,1:i,:)).*squeeze(DPsi_n(1,j,i:-1:1,:))))/i;
        DPsi_n(3,j,i+1,:)=DPsi_n(4,j,i,:)/i;
        DPsi_n(4,j,i+1,:)=sum(squeeze(DPsi_n(2,j,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:))...
            +squeeze(Psi_n(2,1:i,:)).*squeeze(DPsi_n(3,j,i:-1:1,:)))/(2*i);

    end
end
%computing the second partial derivatives on the stadium


D2Psi_n = iv(zeros(4,4,4,max_N+1,length(mu_stadium)));
D2Psi_n(2,1,1,2,:) = inverse_nu_stadium;
D2Psi_n(4,2,3,2,:) = half;
D2Psi_n(4,3,2,2,:) = half;

for i = 2:max_N
    for j = 1:4
        for k = 1:4
            
            D2Psi_n(1,j,k,i+1,:) = D2Psi_n(2,j,k,i,:)/i;
            D2Psi_n(2,j,k,i+1,:) = inverse_nu_stadium.*(-transpose(squeeze( ...
                D2Psi_n(2,j,k,i,:)))+sum(squeeze(D2Psi_n(1,j,k,1:i,:)).*squeeze(Psi_n(1,i:-1:1,:)) + ...
                squeeze(DPsi_n(1,j,1:i,:)).*squeeze(DPsi_n(1,k,i:-1:1,:))))/i;
            D2Psi_n(3,j,k,i+1,:) = D2Psi_n(4,j,k,i,:)/i;
            D2Psi_n(4,j,k,i+1,:) = sum(squeeze(D2Psi_n(2,j,k,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:)) + ...
                squeeze(DPsi_n(2,j,1:i,:)).*squeeze(DPsi_n(3,k,i:-1:1,:)) + ...
                squeeze(DPsi_n(2,k,1:i,:)).*squeeze(DPsi_n(3,j,i:-1:1,:)) + ...
                squeeze(Psi_n(2,1:i,:)).*squeeze(D2Psi_n(3,j,k,i:-1:1,:)))/(2*i);
        end
    end
end

% obtain a uniform bound M on the modulus of the Psi(1:4,1:max_N+1,\mu),
% DPsi(1:4,1:max_N+1,\mu) D2Psi(1:4,1:max_N+1,\mu) when $\mu$ is on the
% stadium.

M_bound = max(max(max(sup(abs(Psi(:,:,:))))));
M_bound = max(M_bound,max(max(max(max(sup(abs(DPsi(:,:,:,:))))))));
M_bound=iv(M_bound);

M_n_bound = max(max(max(sup(abs(Psi_n(:,:,:))))));
M_n_bound = max(M_n_bound,max(max(max(max(sup(abs(DPsi_n(:,:,:,:))))))));
M_n_bound = max(M_n_bound,max(max(max(max(max(sup(abs(D2Psi_n(:,:,:,:,:)))))))));
M_n_bound=iv(M_n_bound);


J = 1;
interpolation_error = iv(1);
interpolation_n_error = iv(1);
while max(sup(interpolation_error),sup(interpolation_n_error)) > err_tol
    J = J+1;
    interpolation_error = M_bound*L_rho/(pie*D_rho*sinh(eta*(J+1)));
    interpolation_n_error = M_n_bound*L_rho/(pie*D_rho*sinh(eta*(J+1)));
end

%
% interpolation nodes
%

theta = (2*iv(0:1:J-1)+1)*pie/(2*J);
mu_tilde = cos(theta);
mu = (mu_left+mu_right)/2+(mu_right-mu_left)*mu_tilde/2;
inverse_nu_grid=mu.^2./(1-mu);
nu_grid = (1-mu)./mu.^2;


% Initialize coefficients of series
Psi = iv(zeros(4,max_N+1,length(mu)));
DPsi = iv(zeros(4,4,max_N+1,length(mu)));

Psi_n = iv(zeros(4,max_N+1,length(mu)));
DPsi_n = iv(zeros(4,4,max_N+1,length(mu)));
D2Psi_n = iv(zeros(4,4,4,max_N+1,length(mu)));

% initial conditions for the profile ODE
Psi(1,1,:) = phi_fun(mu.');
Psi(2,1,:) = w_fun(mu.');
Psi(3,1,:) = u_fun(mu.');
Psi(4,1,:) = z_fun(mu.');

Psi_n(1,1,:) = add_neighborhood(phi_fun(mu),r_NK);
Psi_n(2,1,:) = add_neighborhood(w_fun(mu),r_NK);
Psi_n(3,1,:) = add_neighborhood(u_fun(mu),r_NK);
Psi_n(4,1,:) = add_neighborhood(z_fun(mu),r_NK);

Psi(1,2,:) = Psi(2,1,:);
Psi(2,2,:) = inverse_nu_grid.*transpose(squeeze(-Psi(2,1,:)+half*(Psi(1,1,:).^2-1)));
Psi(3,2,:) = Psi(4,1,:);
Psi(4,2,:) = half*Psi(2,1,:).*Psi(3,1,:);

Psi_n(1,2,:) = Psi_n(2,1,:);
Psi_n(2,2,:) = inverse_nu_grid.*transpose(squeeze(-Psi_n(2,1,:)+half*(Psi_n(1,1,:).^2-1)));
Psi_n(3,2,:) = Psi_n(4,1,:);
Psi_n(4,2,:) = half*Psi_n(2,1,:).*Psi_n(3,1,:);


for i = 2:max_N
    Psi(1,i+1,:) = Psi(2,i,:)/i;
    Psi(2,i+1,:) = inverse_nu_grid.*(-transpose(squeeze(Psi(2,i,:)))...
        +half*sum(squeeze(Psi(1,1:i,:)).*squeeze(Psi(1,i:-1:1,:))))/i;
    Psi(3,i+1,:) = Psi(4,i,:)/i;
    Psi(4,i+1,:) = sum(squeeze(Psi(2,1:i,:)).*squeeze(Psi(3,i:-1:1,:)))/(2*i);

    Psi_n(1,i+1,:) = Psi_n(2,i,:)/i;
    Psi_n(2,i+1,:) = inverse_nu_grid.*(-transpose(squeeze(Psi_n(2,i,:)))...
        +half*sum(squeeze(Psi_n(1,1:i,:)).*squeeze(Psi_n(1,i:-1:1,:))))/i;
    Psi_n(3,i+1,:) = Psi_n(4,i,:)/i;
    Psi_n(4,i+1,:) = sum(squeeze(Psi_n(2,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:)))/(2*i);

end


%initializing the jacobian matrix of Psi with respect to Psi0
DPsi(1,1,1,:)=iv(1);
DPsi(2,2,1,:)=iv(1);
DPsi(3,3,1,:)=iv(1);
DPsi(4,4,1,:)=iv(1);

DPsi_n(1,1,1,:)=iv(1);
DPsi_n(2,2,1,:)=iv(1);
DPsi_n(3,3,1,:)=iv(1);
DPsi_n(4,4,1,:)=iv(1);


DPsi(1,2,2,:) = iv(1);
DPsi(2,1,2,:) = inverse_nu_grid.*transpose(squeeze(Psi(1,1,:)));
DPsi(2,2,2,:) = -inverse_nu_grid; 
DPsi(3,4,2,:) = iv(1);
DPsi(4,2,2,:) = Psi(3,1,:)/2;
DPsi(4,3,2,:) = Psi(2,1,:)/2;


DPsi_n(1,2,2,:) = iv(1);
DPsi_n(2,1,2,:) = inverse_nu_grid.*transpose(squeeze(Psi_n(1,1,:))); % WARNING: inverse_nu_grid needs to be interval
DPsi_n(2,2,2,:) = -inverse_nu_grid; 
DPsi_n(3,4,2,:) = iv(1);
DPsi_n(4,2,2,:) = Psi_n(3,1,:)/2;
DPsi_n(4,3,2,:) = Psi_n(2,1,:)/2;


%computing the jacobian matrix on the interpolation nodes
for i=2:max_N
    for j=1:4
        DPsi(1,j,i+1,:)=DPsi(2,j,i,:)/i;
        DPsi(2,j,i+1,:)=inverse_nu_grid.*(-transpose(squeeze(DPsi(2,j,i,:)))...
            +sum(squeeze(Psi(1,1:i,:)).*squeeze(DPsi(1,j,i:-1:1,:))))/i;
        DPsi(3,j,i+1,:)=DPsi(4,j,i,:)/i;
        DPsi(4,j,i+1,:)=sum(squeeze(DPsi(2,j,1:i,:)).*squeeze(Psi(3,i:-1:1,:))...
            +squeeze(Psi(2,1:i,:)).*squeeze(DPsi(3,j,i:-1:1,:)))/(2*i);
        DPsi_n(1,j,i+1,:)=DPsi_n(2,j,i,:)/i;
        DPsi_n(2,j,i+1,:)=inverse_nu_grid.*(-transpose(squeeze(DPsi_n(2,j,i,:)))... % WARNING! IS THE USE OF invere_nu_grid rigorous here?
            +sum(squeeze(Psi_n(1,1:i,:)).*squeeze(DPsi_n(1,j,i:-1:1,:))))/i;
        DPsi_n(3,j,i+1,:)=DPsi_n(4,j,i,:)/i;
        DPsi_n(4,j,i+1,:)=sum(squeeze(DPsi_n(2,j,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:))...
            +squeeze(Psi_n(2,1:i,:)).*squeeze(DPsi_n(3,j,i:-1:1,:)))/(2*i);
    end
end

%computing the second partial derivatives on the interpolation nodes
D2Psi_n(2,1,1,2,:) = inverse_nu_grid;
D2Psi_n(4,2,3,2,:) = half;
D2Psi_n(4,3,2,2,:) = half;


for i = 2:max_N
    for j = 1:4
        for k = 1:4
            D2Psi_n(1,j,k,i+1,:) = D2Psi_n(2,j,k,i,:)/i;
            D2Psi_n(2,j,k,i+1,:) = inverse_nu_grid.*(-transpose(squeeze( ...
                D2Psi_n(2,j,k,i,:)))+sum(squeeze(D2Psi_n(1,j,k,1:i,:)).*squeeze(Psi_n(1,i:-1:1,:)) + ...
                squeeze(DPsi_n(1,j,1:i,:)).*squeeze(DPsi_n(1,k,i:-1:1,:))))/i;
            D2Psi_n(3,j,k,i+1,:) = D2Psi_n(4,j,k,i,:)/i;
            D2Psi_n(4,j,k,i+1,:) = sum(squeeze(D2Psi_n(2,j,k,1:i,:)).*squeeze(Psi_n(3,i:-1:1,:)) + ...
                squeeze(DPsi_n(2,j,1:i,:)).*squeeze(DPsi_n(3,k,i:-1:1,:)) + ...
                squeeze(DPsi_n(2,k,1:i,:)).*squeeze(DPsi_n(3,j,i:-1:1,:)) + ...
                squeeze(Psi_n(2,1:i,:)).*squeeze(D2Psi_n(3,j,k,i:-1:1,:)))/(2*i);
        end
    end
end


cf_phi = compute_Chebyshev_coeffs(squeeze(Psi(1,:,:)),theta);
abs_phi_bound= sum(abs(cf_phi))+interpolation_error;

mu_int=infsup(inf(mu_left),sup(mu_right));
nu_int=(1-mu_int)./mu_int.^2;
alpha_a=[];
for m=0:1:floor((max_N-2)/2)
    alpha=iv(1);
    for i=1:m+2
        alpha=max(sup(alpha),sup((nu_int^i*abs_phi_bound(i)/2/(m+1))^(1/(i+1))));
    end
    alpha_a=[alpha_a;alpha];
end

[alpha,m]=min(alpha_a);
alpha=iv(alpha);
beta=sqrt(4*alpha + 1)/2 + half;

m=m-1;
C0=iv(0);

cf_u = compute_Chebyshev_coeffs(squeeze(Psi(3,:,:)),theta);
abs_u_bound= sum(abs(cf_u))+interpolation_error;

for i=1:m+1
    C0=iv(max(sup(C0),sup(nu_int^(i-1)*abs_u_bound(i)*(alpha*beta)^(1-i))));
end

N1=2*m+2;

C1_lemma10=iv(zeros(4,1));
for j=1:4
    cf_DPsi1j= compute_Chebyshev_coeffs(squeeze(DPsi(1,j,:,:)),theta);
    abs_cf_DPsi1j= sum(abs(cf_DPsi1j))+interpolation_error;
    for i=1:(N1+1)
        C1_lemma10(j)=iv(max(sup(C1_lemma10(j)),sup(nu_int^(i-1)*alpha^(1-i)*abs_cf_DPsi1j(i))));
    end
end
C2_lemma10=C0*C1_lemma10*nu_int/2/alpha^2;

N2=m+1;
for j=1:4
    for i=1:N2+1
        cf_DPsi3j= compute_Chebyshev_coeffs(squeeze(DPsi(3,j,:,:)),theta);
        abs_cf_DPsi3j= sum(abs(cf_DPsi3j))+interpolation_error;
        C2_lemma10(j)=iv(max(sup(C2_lemma10(j)),sup(nu_int^(i-1)*(alpha*beta)^(1-i)*abs_cf_DPsi3j(i))));
    end
end


%adding error bounds
x1=abs(alpha*delta_x/nu_int);
x2=abs(alpha*beta*delta_x/nu_int);
if sup(x2)>1
    error('x2 is great than 1')
end
N=max_N-1;

best_possible_delta_x = inf(nu_int/(alpha*beta));

d.poly_phi=add_error_bounds(cf_phi,delta_x,interpolation_error,sup(2*(m+1)*alpha^2/nu_int*x1^(N+2)/(1-x1)));
d.poly_w=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi(2,:,:)),theta),delta_x,interpolation_error,sup(2*(m+1)*alpha^3/nu_int^2*x1^(N+2)*(N-2*x1-N*x1+3)/(x1-1)^2));
d.poly_u=add_error_bounds(cf_u,delta_x,interpolation_error,sup(C0*x2^(N+2)/(1-x2)));
d.poly_z=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi(4,:,:)),theta),delta_x,interpolation_error,sup(C0*alpha*beta/nu_int*x2^(N+2)*(N-2*x2-N*x2+3)/(x2-1)^2));
d.poly_DPsi=iv(zeros(4,4,J));

for j=1:4
    d.poly_DPsi(1,j,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(DPsi(1,j,:,:)),theta),delta_x,interpolation_error,sup(C1_lemma10(j)*x1^(N+2)/(1-x1)));
    d.poly_DPsi(2,j,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(DPsi(2,j,:,:)),theta),delta_x,interpolation_error,sup(C1_lemma10(j)*alpha/nu_int*x1^(N+2)*(N-2*x1-N*x1+3)/(x1-1)^2));
    d.poly_DPsi(3,j,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(DPsi(3,j,:,:)),theta),delta_x,interpolation_error,sup(C2_lemma10(j)*x2^(N+2)/(1-x2)));
    d.poly_DPsi(4,j,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(DPsi(4,j,:,:)),theta),delta_x,interpolation_error,sup(C2_lemma10(j)*alpha*beta/nu_int*x2^(N+2)*(N-2*x2-N*x2+3)/(x2-1)^2));
end


d.alpha=alpha;
d.beta=beta;
d.C1_lemma10 = C1_lemma10;
d.C2_lemma10 = C2_lemma10;

d.delta_x = delta_x;
d.interpolation_error = interpolation_error;

cf_phi_n = compute_Chebyshev_coeffs(squeeze(Psi_n(1,:,:)),theta);
abs_phi_bound_n= sum(abs(cf_phi_n))+interpolation_n_error;


alpha_a_n=[];
for m=0:1:floor((max_N-2)/2)
    alpha_n=iv(1);
    for i=1:m+2
        alpha_n=max(sup(alpha_n),sup((nu_int^i*abs_phi_bound_n(i)/2/(m+1))^(1/(i+1))));
    end
    alpha_a_n=[alpha_a_n;alpha_n];
end

[alpha_n,m]=min(alpha_a_n);
alpha_n=iv(alpha_n);
beta_n=sqrt(4*alpha_n + 1)/2 + half;

m=m-1;
C0_n=iv(0);

cf_u_n = compute_Chebyshev_coeffs(squeeze(Psi_n(3,:,:)),theta);
abs_u_n_bound= sum(abs(cf_u_n))+interpolation_n_error;

for i=1:m+1
    C0_n=iv(max(sup(C0_n),sup(nu_int^(i-1)*abs_u_n_bound(i)*(alpha_n*beta_n)^(1-i))));
end

N1=2*m+2;

C1_lemma10_n=iv(zeros(4,1));
for j=1:4
    cf_DPsi1j_n= compute_Chebyshev_coeffs(squeeze(DPsi_n(1,j,:,:)),theta);
    abs_cf_DPsi1j_n= sum(abs(cf_DPsi1j_n))+interpolation_n_error;
    for i=1:(N1+1)
        C1_lemma10_n(j)=iv(max(sup(C1_lemma10_n(j)),sup(nu_int^(i-1)*alpha_n^(1-i)*abs_cf_DPsi1j_n(i))));
    end
end
C2_lemma10_n=C0_n*C1_lemma10_n*nu_int/2/alpha_n^2;

N2=m+1;
for j=1:4
    for i=1:N2+1
        cf_DPsi3j_n= compute_Chebyshev_coeffs(squeeze(DPsi_n(3,j,:,:)),theta);
        abs_cf_DPsi3j_n= sum(abs(cf_DPsi3j_n))+interpolation_n_error;
        C2_lemma10_n(j)=iv(max(sup(C2_lemma10_n(j)),sup(nu_int^(i-1)*(alpha_n*beta_n)^(1-i)*abs_cf_DPsi3j_n(i))));
    end
end
N3=2*m+3;
C1_lemma11_n=iv(zeros(4,4));
C2_lemma11_n=iv(zeros(4,4));
for i=1:4
    for j=1:4
        C1_lemma11_n(i,j)=C1_lemma10_n(i)*C1_lemma10_n(j)*nu_int/alpha_n/(2*alpha_n-1);
        cf_D2Psi1ij_n= compute_Chebyshev_coeffs(squeeze(D2Psi_n(1,i,j,:,:)),theta);
        abs_cf_D2Psi1ij_n= sum(abs(cf_D2Psi1ij_n))+interpolation_n_error;
        for k=1:(N3+1)
            C1_lemma11_n(i,j)=iv(max(sup(C1_lemma11_n(i,j)),sup(nu_int^(k-1)*alpha_n^(1-k)*abs_cf_D2Psi1ij_n(k))));
        end
    end
end



for i=1:4
    for j=1:4
        C2_lemma11_n(i,j)=(C1_lemma10_n(i)*C2_lemma10_n(j)+C1_lemma10_n(j)*C2_lemma10_n(i)+C0_n*C1_lemma11_n(i,j))/(2*alpha_n^2);
        cf_D2Psi3ij_n= compute_Chebyshev_coeffs(squeeze(D2Psi_n(3,i,j,:,:)),theta);
        abs_cf_D2Psi3ij_n= sum(abs(cf_D2Psi3ij_n))+interpolation_n_error;
        for k=1:N2+1
            C2_lemma11_n(i,j)=iv(max(sup(C2_lemma11_n(i,j)),sup(nu_int^(k-1)*(alpha_n*beta_n)^(1-k)*abs_cf_D2Psi3ij_n(k))));
        end
    end
end


%adding error bounds
x1_n=abs(alpha_n*delta_x/nu_int);
x2_n=abs(alpha_n*beta_n*delta_x/nu_int);
if sup(x2_n)>1
    error('x2_n is great than 1')
end
N=max_N-1;
NN_j=101;
delta_x_a = linspace(0,sup(delta_x),NN_j);
delta_x_a = iv(delta_x_a(1:end-1),delta_x_a(2:end));
mu_a= linspace(inf(mu_left),sup(mu_right),NN);
mu_a = iv(mu_a(1:end-1),mu_a(2:end));
u_n=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi_n(3,:,:)),theta),delta_x_a(1),interpolation_n_error,sup(C0_n*x2_n^(N+2)/(1-x2_n)));
u_n_fun = @(mu)eval_cf_z(mu,u_n,mu_left,mu_right);
%check if u is sign-definite for mu in [mu_left,mu_right] on [delta_x 0]
u_n_i=u_n_fun(mu_a(1));
if in(0,u_n_i)
    d.u_dx_vs_0=0;
elseif inf(u_n_i)>0
    d.u_dx_vs_0=1;
else
    d.u_dx_vs_0=-1;
end
i=2;
while d.u_dx_vs_0~=0 && i<NN
    u_n_i=u_n_fun(mu_a(i));
    if in(0,u_n_i)
        d.u_dx_vs_0=0;
    elseif (inf(u_n_i)>0)*d.u_dx_vs_0==-1
        d.u_dx_vs_0=0;
    end
    i=i+1;
end
j=2;
while d.u_dx_vs_0~=0 && j<NN_j
u_n=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi_n(3,:,:)),theta),delta_x_a(j),interpolation_n_error,sup(C0_n*x2_n^(N+2)/(1-x2_n)));
u_n_fun = @(mu)eval_cf_z(mu,u_n,mu_left,mu_right);
i=1;
while d.u_dx_vs_0~=0 && i<NN
    u_n_i=u_n_fun(mu_a(i));
    if in(0,u_n_i)
        d.u_dx_vs_0=0;
    elseif (inf(u_n_i)>0)*d.u_dx_vs_0==-1
        d.u_dx_vs_0=0;
    end
    i=i+1;
end
j=j+1;
end
z_n=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi_n(4,:,:)),theta),delta_x_a(1),interpolation_n_error,sup(C0_n*alpha_n*beta_n/nu_int*x2_n^(N+2)*(N-2*x2_n-N*x2_n+3)/(x2_n-1)^2));
z_n_fun = @(mu)eval_cf_z(mu,z_n,mu_left,mu_right);
z_n_i=z_n_fun(mu_a(1));
if in(0,z_n_i)
    d.z_dx_vs_0=0;
elseif inf(z_n_i)>0
    d.z_dx_vs_0=1;
else
    d.z_dx_vs_0=-1;
end
i=2;
while d.z_dx_vs_0~=0 && i<NN
    z_n_i=z_n_fun(mu_a(i));
    if in(0,z_n_i)
        d.z_dx_vs_0=0;
    elseif (inf(z_n_i)>0)*d.z_dx_vs_0==-1
        d.z_dx_vs_0=0;
    end
    i=i+1;
end
j=2;
while d.z_dx_vs_0~=0 && j<NN_j
z_n=add_error_bounds(compute_Chebyshev_coeffs(squeeze(Psi_n(4,:,:)),theta),delta_x_a(j),interpolation_n_error,sup(C0_n*alpha_n*beta_n/nu_int*x2_n^(N+2)*(N-2*x2_n-N*x2_n+3)/(x2_n-1)^2));
z_n_fun = @(mu)eval_cf_z(mu,z_n,mu_left,mu_right);
i=1;
while d.z_dx_vs_0~=0 && i<NN
    z_n_i=z_n_fun(mu_a(i));
    if in(0,z_n_i)
        d.z_dx_vs_0=0;
    elseif (inf(z_n_i)>0)*d.z_dx_vs_0==-1
        d.z_dx_vs_0=0;
    end
    i=i+1;
end
j=j+1;
end


for j=1:4
    for i=1:4
        d.poly_D2Psi_n(1,j,i,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(D2Psi_n(1,j,i,:,:)),theta),delta_x,interpolation_n_error,sup(C1_lemma11_n(j,i)*x1_n^(N+2)/(1-x1_n)));
        d.poly_D2Psi_n(2,j,i,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(D2Psi_n(2,j,i,:,:)),theta),delta_x,interpolation_n_error,sup(C1_lemma11_n(j,i)*alpha_n/nu_int*x1_n^(N+2)*(N-2*x1_n-N*x1_n+3)/(x1_n-1)^2));
        d.poly_D2Psi_n(3,j,i,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(D2Psi_n(3,j,i,:,:)),theta),delta_x,interpolation_n_error,sup(C2_lemma11_n(j,i)*x2_n^(N+2)/(1-x2_n)));
        d.poly_D2Psi_n(4,j,i,:)=add_error_bounds(compute_Chebyshev_coeffs(squeeze(D2Psi_n(4,j,i,:,:)),theta),delta_x,interpolation_n_error,sup(C2_lemma11_n(j,i)*alpha_n*beta_n/nu_int*x2_n^(N+2)*(N-2*x2_n-N*x2_n+3)/(x2_n-1)^2));
    end
end


d.alpha_n=alpha_n;
d.beta_n=beta_n;


d.C1_lemma11_n = C1_lemma11_n;
d.C2_lemma11_n = C2_lemma11_n;
d.delta_x = delta_x;
d.interpolation_n_error = interpolation_n_error;




d.run_time = toc(start_time);


