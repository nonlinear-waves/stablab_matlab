function out = bc_fun_cheby(fun,L,n)

ya = fun(-L);
yb = fun(L);


PL = [-sqrt(1+3*n),1];
PR = [sqrt(1+3*n),1];

out = [PL*ya;PR*yb];











