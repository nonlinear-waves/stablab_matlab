function out = A_coord_change(x,lambda,s,p)
% 
% Eigenvalue problem: lambda*u + (Au)' = u'' in integrated coordinates

% bar_v = -p.gamma*tanh(p.gamma*x/2);
 
bar_v_prime = -0.5*p.gamma^2*(1-tanh(p.gamma*x/2)^2);

q = p.M*pi*bar_v_prime;

out = [ 0 q 1 0; -q 0 0 1; lambda 0 1 q; 0 lambda -q -1];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out = rot(theta)
% 
% out = [cos(theta) -sin(theta); sin(theta) cos(theta)];

