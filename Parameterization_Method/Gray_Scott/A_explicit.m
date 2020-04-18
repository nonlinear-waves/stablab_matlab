function out = A_explicit(x,lambda,s,p)

% Evans function matrix
%
% W'(x) = A(x,\lambda)W(x)
%
% W = (u,u',v,v')^T

% only valid if gamma*alpha = 1 and 0 < gamma < 2/9
u = 1-3*p.gamma./(1+p.Q*cosh(x/sqrt(p.gamma)));
v =3./(1+p.Q*cosh(x/sqrt(p.gamma)));

out = [0 1 0 0; 
    lambda+v^2+p.alpha 0 2*u*v 0;
    0 0 0 1;
    -v^2/p.gamma 0 lambda+(1-2*u*v)/p.gamma 0];



