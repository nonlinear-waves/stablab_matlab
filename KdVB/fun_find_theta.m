function out = fun_find_theta(u,fun,y0)

% projection onto the first two components minus y0.

temp = fun(u(1),u(2));

out = temp(1:2)-y0;