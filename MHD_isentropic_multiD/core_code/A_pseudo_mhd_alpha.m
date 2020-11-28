function out = A_pseudo_mhd_alpha(x,lambda,s,p)

if strcmp(s.sol.solver,'bvp_fsolve')
    temp = s.sol.deval(x);
else
    temp = soln(x,s);
end

out = temp(1)*A_multi2(x,lambda,s,p);







