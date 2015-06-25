function out=capa(x,y,lambda,s,p,A,n,k,MU)
% out=capa(x,y,lambda,s,p,A,n,k,MU)
%
% Returns the value y'(x) of the first order system y'=(A(x,lambda)-(mu)*Identity)*y
%
% Input "x" is the value where y'(x) is evaluated, $y$ is the vector y(x),
% "lambda" is the point in the complex plane where the Evans function is
% evaluated, s,p are structures explained in the STABLAB documentation, "A"
% is a function handle to the Evans matrix, "n" is the dimension of the 
% system and "k" is the dimension of the manifold sought, and "MU" is the
% eigenvalue corresponding to the largest or smallest eigenvalue of
% A(\pm \infty,lambda) 

out = (A(x,lambda,s,p)-MU*eye(nchoosek(n,k))) * y;
