function out = bc(ya,yb,n,lam,p)

% We get eigenvalues for A+ = A-,
% lambda = -
%
%
%
mu1 = sqrt(lam*n+p.alpha);
mu2 = sqrt(lam*n+1/p.gamma);

PL = [mu1,-1,0,0;
      0, 0, mu2, -1];
  
PR = [-mu1,-1,0,0;
      0, 0, -mu2,-1];




out = [PL*ya;PR*yb];











