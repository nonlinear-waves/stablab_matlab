function out = A(x,lambda,s,p)

v = (soln(x,s));

ee=-0;

out=[0+ee, 1, 0;
lambda-v(2)/v(1)^2*exp(-1/v(1)), -p.c+ee, -exp(-1/v(1));
p.be/p.c*v(2)/v(1)^2*exp(-1/v(1)),0,lambda/p.c+p.be/p.c*exp(-1/v(1))+ee];






