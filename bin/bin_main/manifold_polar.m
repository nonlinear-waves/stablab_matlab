function [Omega,gamma] = manifold_polar(x,y,lambda,A,s,p,m,k,mu)
% [omega,gamma] = manifold_polar(x,y,lambda,A,s,p,m,k)
%
% Returns "Omega", the orthogonal basis for the manifold evaluated at x(2)
% and "gamma" the radial equation evaluated at x(2).
%
% Input "x" is the interval on which the manifold is solved, "y" is the
% initializing vector, "lambda" is the point in the complex plane where the
% Evans function is evaluated, "A" is a function handle to the Evans
% matrix, s, p,and m are structures explained in the STABLAB documentation, 
% and k is the dimension of the manifold sought.

% Solve the ODE
[unused,Y] = m.ode_fun(m.method,[x(1) x(2)],[reshape(y,m.n*k,1);0],m.options, ...
    lambda,A,s,p,m.n,k,mu);

% Reshape the (n*k+1) vector to be the n x k matrix Omega
Omega = reshape(Y(end,1:m.n*k).',m.n,k);

% gamma = e^{rho}, rho = log(gamma)
gamma = exp(Y(end,m.n*k+1));


