function out = A_lop_alpha(x,lambda,s,p)

if x<0
    rho = p.rho_m;
    u1 = p.u1_m;
    p_rho = p.a0*p.gamma*p.rho_m^(p.gamma-1);
elseif x >= 0
    rho = p.rho_p;
    u1 = p.u1_p;
    p_rho = p.a0*p.gamma*p.rho_p^(p.gamma-1);
end
h1 = p.h1;

% alpha model with alpha = 1
A0 = [ 1 0 0 0 0;
    u1 rho 0 0 0;
    0 0 rho 0 0;
    0 0 0 1 0;
    0 0 0 0 1];

A1 = [
  u1 rho 0 0 0;
  u1^2+p_rho 2 0 -h1 0;
   0 0 1 0 -h1;
   0 0 0 1 0;
   0 0 -h1 0 u1
];

A2 = [
  0 0 rho 0 0;
  0 0 1 0 -h1;
  p_rho 0 0 h1 0;
  0 0 h1 0 1-u1;
  0 0 0 0 0
];

out = -A1\(lambda*A0+1i*s.xi*A2);











