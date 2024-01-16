function out = bc_u_ode(ya,yb,data)

% boundary conditions for solving the BVP for the profile and Ricatti
% equation.

out = [ya(1:2)-[data.phi0;data.w0];yb(3:4)-[data.u1;data.z1]];