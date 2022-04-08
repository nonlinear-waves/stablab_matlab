function out = J(U,p)

out=[-p.c -p.c/p.be;
     1/U(1)^2*p.be/p.c*U(2)*exp(-1/U(1)) p.be/p.c*exp(-1/U(1))];

