function out = A(x,lambda,s,p)

% profile solution
v = (soln(x,s));

%Integrated system
out=[0, 1, 0, 0, 0, 0; 
    0, 0, 1, 0, 0, 0; 
    0, 0, 0, 1, 0, 0; 
    -3*lambda/(v(1)^3*p.C), (-3*v(4)^2*v(1)^3*p.ss-15*p.D*v(4)*v(1)^2*p.ss-9*v(4)^2*p.K1*v(1)^2-24*p.D^2*v(1)*p.ss-36*p.D*v(4)*p.K1*v(1)-36*p.D^2*p.K1)/(p.C*v(1)^4*(v(4)*v(1)+4*p.D)*(v(4)*v(1)+p.D)),...
    (1/4)*p.be*(v(4)*v(1)+4*p.D)/(p.C*(v(4)*v(1)+p.D)), 0, (3/2)*lambda/(v(1)*p.C*(v(4)*v(1)+p.D)), 3*p.D*(v(1)*p.ss+3*p.K1)/(p.C*v(1)^2*(v(4)*v(1)+4*p.D)*(v(4)*v(1)+p.D)); 
    0, 0, 0, 0, 0, 1; 
    0, 2*v(4)*(2*v(4)*v(1)^2*p.ss+6*p.D*v(1)*p.ss+3*v(4)*p.K1*v(1)+6*p.D*p.K1)/(v(1)^2*(v(4)*v(1)+4*p.D)*(v(4)*v(1)+p.D)), -p.be*v(4)*v(1)^2/(2*v(4)*v(1)+2*p.D), 0,...
    lambda/(v(4)*v(1)+p.D), 2*p.D*(v(1)*p.ss+3*p.K1)/(v(1)*(v(4)*v(1)+4*p.D)*(v(4)*v(1)+p.D))];
