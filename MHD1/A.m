function out=A(x,lambda,s,p)
%MHD parallel case

v=soln(x,s);

out=[...
    0 1/p.mu 0 0;
    lambda*v v/p.mu 0 -p.sigma*p.B*v;
    0 0 0 v*p.sigma*p.mu0;
    0 -p.B*v/p.mu lambda*v v^2*p.sigma*p.mu0...
    ];



