function [rho,rho_x,u1,u1_x,u2,u2_x,h1,h1_x,h2,h2_x] = eval_prof_mhd_alpha(s,p,x)

if strcmp(s.sol.solver,'bvp_fsolve')
    temp = s.sol.deval(x);
else
    temp = soln(x,s);
end


u1 = temp(1);
u1_x = profile_ode_pseudo(x,u1,s,p);

rho = 1/u1;
rho_x = -u1_x/u1^2;

h1 = p.h1;
h1_x = 0;

u2 = 0;
u2_x = 0;

h2 = 0;
h2_x = 0;

    







