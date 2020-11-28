function [s,p] = profile_solve_pseudo2(s,p,num_inf)
% [s,p] = profile(s,p)
%
% solve the profile for isentropic in Eulerian coordinates with
% wave speed zero.

%
% parameters
%

% dependent parameters
p.u1_m = 1;
p.h2_m = 0;
p.h2_p = 0;
p.u2_m = 0;
p.u2_p = 0;
p.rho_m = 1/p.u1_m;
p.rho_p = 1/p.u1_p;
p.alpha = 1;
p.spd = 0;
p.m1 = 1;

p.a = p.u1_p^p.gamma*(1-p.u1_p)/(1-p.u1_p^p.gamma);

p.a0 = p.a;
p.u_plus = p.u1_p;
p.rho_plus = 1/p.u1_p;
p.rho_minus = 1;
p.u_minus = 1;


s.F = @(x,y)(profile_ode_pseudo(x,y,s,p));

s.Flinear = @profile_jacobian;
s.UL = p.u_minus;
s.UR = p.u_plus;

if nargin < 3
s.I = 30;
else
   s.I = num_inf;
end
s.L = -s.I;
s.R = s.I;
s.side = 1;

% solve the ODE
a = 0.5*(s.UL+s.UR);
c = 0.5*(s.UL-s.UR);
solinit.guess = @(x)(a-c*tanh(x));
stats = 'off';
options = bvp_fsolve_options('algorithm_stats',stats);
solinit.boundaries = [s.L,s.R];

bc_fun =@(fun)(fun(0)-a);
sol = bvp_fsolve(s.F,bc_fun,solinit,options);
% err check the solution
[ode_err,bc_err ] = sol.check_err(1001);
sol.ode_err = ode_err;
sol.bc_err = bc_err;

s.sol = sol;




