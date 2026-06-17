function out = A(x,lambda,~,p)

rho_0 = p.rho_m +p.rho_0_der*(x+p.h);

fun = p.rho_0_der/rho_0;

out = [0, 1; p.k^2*(1-p.g*fun/lambda^2), -fun]; 



