function out = A(x,lambda,s,p)

u0 = sqrt(2)*sech(x);
v0 = 0;

% out = [ 0 0 1 0;
%     0 0 0 1;
%     1-3*q^2, lambda+2*p.nu, 0, 0;
%     -lambda, p.mu-q^2, 0, 0];

out = [0,1,0,0;
    1-3*u0^2-v0^2, 0, lambda+2*p.nu-2*u0*v0, 0;
    0, 0, 0, 1;
    -lambda-2*u0*v0, 0, p.mu-3*v0^2-u0^2, 0];
    


