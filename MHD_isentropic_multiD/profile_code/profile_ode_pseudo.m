function out = profile_ode_pseudo(x,u,s,p)
%  out = profile_ode(x,u,s,p)
%
% isentropic profile ODE in Euler coordinates

rho = 1/u;
out = u-p.u_minus+ P(rho,p)-P(p.rho_minus,p);
out = (1/(2*p.mu+p.eta))*out;

% pseudo Lagrangian coordinates
out = u*out;


