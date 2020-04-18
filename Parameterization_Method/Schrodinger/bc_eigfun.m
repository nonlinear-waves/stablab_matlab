function out = bc_eigfun(ya,yb,lambda,p)



Q = 0.5*sqrt((p.mu+1)^2-4*(lambda^2+2*lambda*p.nu+p.mu));
C = 0.5*(1+p.mu);
mu1 = (sqrt(Q+C));
mu2 = (sqrt(C-Q));

PL = [mu1,-1,0,0;
      0, 0, mu2, -1];
  
PR = [-mu1,-1,0,0;
      0, 0, -mu2,-1];
  
out = [PL*ya(1:4);PR*yb(5:8);norm(yb(1:4))^2-1;yb(1:4)-ya(5:8)];



