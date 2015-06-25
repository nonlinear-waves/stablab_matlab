function out = A(x,lambda,s,p)
% out = A(x,lambda,s,p)

%{
<bl>
We write the eigenvalue problem as the first order system,
\eq{
W' = A(x,\lambda)W,
}{}
where
\eqn{
A(x,\lambda) = \mat{0&\lambda&1\\
0&0&1\\
\lambda v&\lambda v& f(v,\gamma,a)-\lambda},
}{}
and
\eqn{
f(v,\gamma,a) = v-v^{-\gamma}(-v^{\gamma+1}+a(\gamma-1)+(a+1)v^{\gamma}).
}{}
<el>
%}

% evaluate profile
v = soln(x,s);

f = v-v^(-p.gamma)*(-v^(p.gamma+1)+p.a*(p.gamma-1)+(p.a+1)*v^(p.gamma));

out = [0                       lambda         1;
          0                       0                   1;
          lambda*v          lambda*v       f-lambda];
      
      
      
      
      
      
  