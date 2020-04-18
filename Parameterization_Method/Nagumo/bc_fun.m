function out = bc_fun(ya,yb,n)


PL = [-sqrt(1+3*n),1];
PR = [sqrt(1+3*n),1];

out = [PL*ya;PR*yb];











