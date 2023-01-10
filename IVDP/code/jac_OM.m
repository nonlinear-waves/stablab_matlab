function out = jac_OM(v,pars)
% assumes that f_xy = f_yx and g_xy = g_yx

x = v(1);
y = v(2);
p = v(3);
q = v(4);

d = fg(x,y,pars);

out = [
    d.f_x, d.f_y, pars.sigma1^2, 0;
    d.g_x, d.g_y, 0, pars.sigma2^2;
    -p*d.f_xx-q*d.g_xx, -p*d.f_xy-q*d.g_xy, -d.f_x, -d.g_x;
    -p*d.f_xy-q*d.g_xy, -p*d.f_yy-q*d.g_yy, -d.f_y, -d.g_y
];

out(3,1) = out(3,1)+4*pars.a*pars.eps;