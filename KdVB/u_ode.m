function out = u_ode(x,y,nu)

% the RHS of the ODE for the BVP. Components are (phi,w,u,z).
out = [y(2);(1/nu)*(-y(2)+0.5*(y(1)^2-1))  ;y(4); y(2)*y(3)/2];