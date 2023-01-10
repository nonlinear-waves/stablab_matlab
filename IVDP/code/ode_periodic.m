function out = ode_periodic(~,v,X,pars)

x = v(1);
y = v(2);

d = fg(x,y,pars);

out = X*[d.f; d.g];











