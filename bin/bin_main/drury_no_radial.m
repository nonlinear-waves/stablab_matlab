function ydot=drury_no_radial(t,y,lambda,A,s,p,n,k,not_used_1,not_used_2)
% ydot=drury(t,y,lambda,A,s,p,n,k,mu,damping)
%
% Returns the ODE output for the polar method using the method of Drury
%
% Input "t" and "y" are provided by ode45, "A" is a function handle to the
% desired Evans matrix, s,p are structures explained in the STABLAB
% documentation, "n" is the dimension of the system and "k" is the
% dimension of the manifold.

W = reshape(y(1:k*n,1),n,k);
A_temp = A(t,lambda,s,p);
ydot = [reshape((eye(n)-W*W')*A_temp*W,n*k,1);0];
