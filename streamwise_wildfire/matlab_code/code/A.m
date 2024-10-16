function out = A(x,lambda,s,p)

temp = soln(x,s);

u = temp(1);

[~,df] = f(u,p);


out = [0, 1; lambda-df, Wz(x,p)-p.c];
