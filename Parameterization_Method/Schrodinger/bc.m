function out = bc(ya,yb,n,m,lam1,lam2,p)

psi1 = sqrt(0.5*(1+p.mu+sqrt((1+p.mu)^2-4*(p.mu+(lam1*n+lam2*m)*(lam1*n+lam2*m+2*p.nu)))));
psi2 = sqrt(0.5*(1+p.mu-sqrt((1+p.mu)^2-4*(p.mu+(lam1*n+lam2*m)*(lam1*n+lam2*m+2*p.nu)))));

psi1_p = psi1;
psi1_m = -psi1;
psi2_p = psi2;
psi2_m = -psi2;

A_pm = [0 1 0 0;
    1 0 (lam1*n+lam2*m)+2*p.nu 0;
    0 0 0 1;
    -(lam1*n+lam2*m) 0 p.mu 0];
[V,D] = eig(A_pm');

vec = abs(diag(D)-conj(psi1_p));
ind1 = find(vec == min(vec));
vec = abs(diag(D)-conj(psi2_p));
ind2 = find(vec == min(vec));

vec = abs(diag(D)-conj(psi1_m));
ind3 = find(vec == min(vec));
vec = abs(diag(D)-conj(psi2_m));
ind4 = find(vec == min(vec));

% if (n==1)&&(m==1)
%    
%     PR = V(:,[ind1,ind4])';
%     PL = V(:,[ind2,ind3])'; 
%     
% else
    PR = V(:,[ind1,ind2])';
    PL = V(:,[ind3,ind4])';
% end

out = [PL*ya;PR*yb];

