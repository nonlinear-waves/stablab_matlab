function s = solve_profile(s,p,speed_guess,guess_fun)


s.L = -s.I;
s.R = s.I;
s.side = 1;
s.rarray = [1,2];
s.larray = [3,4];
s.phase = p.alpha*s.UL(1)+(1-p.alpha)*s.UR(1);
s.bvp_options = bvpset('RelTol', 1e-10, 'AbsTol', 1e-10,'Nmax', 20000);



ode = @(x,y,c) [profile_ode(x,y(s.rarray),c,p); ...
    -profile_ode(-x,y(s.larray),c,p)];

bc = @(ya,yb,c)bc_fun(ya,yb,c,p,s);

x_dom = linspace(0,s.R,30);



solinit = bvpinit(x_dom,guess_fun,speed_guess);

s.sol = bvp5c(ode,bc,solinit,s.bvp_options);