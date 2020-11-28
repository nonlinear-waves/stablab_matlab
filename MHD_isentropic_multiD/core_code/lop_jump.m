function [f0,ftilde,A1] = lop_jump(x,s,p)

% variables
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

pr = p.a0*rho^p.gamma;

%f0
f0 = [ ...
rho;
1;
0;
h1;
0
];

%ftilde
ftilde = s.xi*[
    0;
    0;
    pr+0.5*h1^2;
    0;
    0
];

% alpha = 1
A1 = [
  u1 rho 0 0 0;
  u1^2+p_rho 2 0 -h1 0;
   0 0 1 0 -h1;
   0 0 0 1 0;
   0 0 -h1 0 u1
];






