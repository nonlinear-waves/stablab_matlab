function out = A(x,lambda,s,p)

v = soln(x,s);
b = 1/v(1);
h = 1-p.a*p.gamma*v(1)^(-p.gamma-1) + v(2)/v(1)^2-lambda/v(1);

out = [0        lambda   1        0;
       0        0        1        0;
       0        0        0        1;
       lambda/p.d lambda/p.d h/p.d      (-b/p.d)];

