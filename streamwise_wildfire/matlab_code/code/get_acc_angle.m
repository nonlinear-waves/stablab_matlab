function out = get_acc_angle(omega,s,p)

ode = @(x,theta) angle_ode(x,theta,omega,s,p);

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% eigenvalues of A = [0, 1; -dfv+omega,-p.c+Wz(-infty,p)]
[~,dfv] = f(s.UL(1),p);
eg_LH_pos = (Wz(-1,p)-p.c+sqrt((p.c-Wz(-1,p))^2-4*(dfv-omega)))/2;
% eg_LH_neg = (Wz(-1,p)-p.c-sqrt((p.c-Wz(-1,p))^2-4*(dfv-omega)))/2;

% eigenvalues of A = [0, 1; -dfv+omega,-p.c+Wz(+infty,p)]
[fv,dfv] = f(s.UR(1),p);
% eg_RH_pos = (Wz(1,p)-p.c+sqrt((p.c-Wz(1,p))^2-4*(dfv-omega)))/2;
eg_RH_neg = (Wz(1,p)-p.c-sqrt((p.c-Wz(1,p))^2-4*(dfv-omega)))/2;

sol_angle = ode15s(ode,[s.L,s.R],atan(eg_LH_pos),options);

out = sol_angle.y(end);

% out = mod(sol_angle.y(end)-eg_RH_neg-pi,2*pi)-pi;
