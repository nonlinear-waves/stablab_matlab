function out = Ak(x,lambda,s,p)

gamma = .5*sqrt(1-p.S^2);
u = 1.5*(1-p.S^2)*(sech(gamma*x)).^2;
ux = -2*gamma*u.*tanh(gamma*x);
uxx = 2*gamma^2*u.*(2-3*(sech(gamma*x)).^2);

a41 = -lambda^2-2*uxx;
a42 = 2*lambda*p.S-4*ux;
a43 = (1-p.S^2)-2*u;

out =  [0      1      0      0      0       0; 
        0      0      1      1      0       0;
        a42    a43    0      0      1       0;
        0      0      0      0      1       0;
        (-a41) 0      0      a43    0       1;
        0      (-a41) 0      (-a42) 0       0];
