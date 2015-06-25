function out = J(U,p)


out = [0 1;(1+p.a*p_prime(U(1),p))/p.d 1/(-U(1)*p.d)];

% -----------------------------------------------------------
function out = p_prime(v,p)


out = -p.gamma * v^(-p.gamma-1);