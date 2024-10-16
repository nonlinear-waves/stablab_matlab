function out = angle_ode(x,theta,omega,s,p)

temp = soln(x,s);
u = temp(1); 

[~,df] = f(u,p);

out = (omega-df)*cos(theta)^2+(Wz(x,p)-p.c)*cos(theta)*sin(theta)-sin(theta)^2;

