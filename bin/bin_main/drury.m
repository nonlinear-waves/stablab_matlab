function ydot = drury(t,y,lambda,A,s,p,n,k,mu)
% ydot = drury(t,y,lambda,A,s,p,n,k,mu,damping)
%
% Returns the ODE output for the polar method using the method of Drury
%
% Input "t" and "y" are provided by the ode soler, "A" is a function handle to the
% desired Evans matrix, s and p are structures explained in the STABLAB
% documentation, "n" is the dimension of the system and "k" is the
% dimension of the manifold

% W = Omega*alpha, gamma = det(alpha), rho = log(gamma)

% Reshape the (k*n+1) vector to be the n x k matrix Omega 
Omega = reshape(y(1:k*n,1),n,k);

% Evaluate A(x,lambda)
A_temp = A(t,lambda,s,p);

% Compute Omega' and rho'
ydot = [reshape((eye(n)-Omega*Omega')*A_temp*Omega,n*k,1); ...
    trace(Omega'*A_temp*Omega)-mu];

