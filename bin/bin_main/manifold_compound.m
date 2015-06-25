function out = manifold_compound(x,z,lambda,s,p,m,A,k,pmMU)
% out = manifold_compound(x,z,lambda,s,p,m,A,k,pmMU)
%
% Returns the vector representing the manifold evaluated at x(2).
%
% Input "x" is the interval the manifold is computed on, "z" is the
% initializing vector for the ode solver, "lambda" is the point on the
% complex plane where the Evans function is computed, s,p,m are structures
% explained in the STABLAB documentation, "A" is the function handle to the
% desired Evans matrix, "k" is the dimension of the manifold sought, and
% "pmMU" is 1 or -1 depending on if respectively the growth or decay 
% manifold is sought.

[R,D] = eig(A(x(1),lambda,s,p));
[e,mat] = max(real(diag(pmMU*D)));
MU = D(mat,mat);

[X,Z]=ode45(@capa,x,z,m.options,lambda,s,p,A,m.n,k,MU);
out = Z(end,:).';