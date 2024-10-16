function RHS = profile_ode(x,y,c,p)

RHS = [y(2);-f(y(1),p)-(c-Wz(x,p))*y(2)];

