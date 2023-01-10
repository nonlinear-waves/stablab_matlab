function d = fg(x,y,p)

a = p.a;

d.f = y;
d.g = 2*a*y*(x^2 - 1) - x;
d.f_x = 0;
d.f_y = 1;
d.g_x = 4*a*x*y - 1;
d.g_y = 2*a*(x^2 - 1);
d.f_xx = 0;
d.f_xy = 0;
d.f_yy = 0;
d.g_xx = 4*a*y;
d.g_xy = 4*a*x;
d.g_yy = 0;




