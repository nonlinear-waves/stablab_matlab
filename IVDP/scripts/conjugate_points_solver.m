clc; close all; clear all;

%
% system parameters
%

p.a = 0.5;
p.eps = 0;
p.sigma1 = 1; % fixed to be 1
p.sigma2 = 1; % fixed to be 1
p.eps = 0; % eps = 0 corresponds to FW.

theta = 3; % manifold angle
r_theta = 1; % radius of circle input
L = 35;

%
% load the profile and the unstable manifold
%

curr_dir = cd; cd('../data');
ld = load(replace(['manifold_a_',num2str(p.a),'_eps_',num2str(p.eps)],'.','P'));
UM = ld.UM;
ld = load(replace(['profile_a_',num2str(p.a)],'.','P'));
sol_per = ld.sol_per;
cd('../code');

% 
% initialize the Plucker coordinates
%

A = jac_OM(zeros(4,1),p);
[V,D] = eigs(A);
egs = diag(D);
ind = find(diag(D)>0);
mu1 = D(ind(1),ind(1));
mu2 = D(ind(2),ind(2));
muI = (mu1+mu2)*eye(6);
Vbasis = V(:,ind);
y0E0 = wedgie([real(Vbasis(:,1)),imag(Vbasis(:,2))]);

theta1 = r_theta*(cos(theta)+1i*sin(theta));
theta2 = r_theta*(cos(theta)-1i*sin(theta));

y0 = real(eval_mani(UM,theta1,theta2)); % point on the manifold circle with angle theta

% go back along the path to get close to the fixed point
y0_back = real(eval_mani(UM,theta1*exp(-mu1*L),exp(-mu2*L)*theta2)); 

options = odeset('RelTol',1e-13,'AbsTol',1e-20);

ode_fun = @(t,y) ode_hamiltonian_bound_OM(t,y,p);

L = 35;

% solve for river path forward
lastwarnd = '';
sol_mani = ode15s(ode_fun,[0,20],y0,options);

% solve for river path backward piece
lastwarn = '';
sol_mani_backward = ode15s(ode_fun,[0,L],y0_back,options);

% do error checking 
err = norm(deval(sol_mani_backward,L)-y0)/norm(y0);
if err > 1e-5
    err
    error('error detected');
end

% Find any conjugate points
[conj_forward,conj_backward] = get_conjugate_points(p,sol_mani,sol_mani_backward,muI,y0E0);

cd('../scripts');