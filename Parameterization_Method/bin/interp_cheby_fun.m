function out = interp_cheby_fun(cf,x,N,a,b)


x0 = (x - 0.5*(a+b))/(0.5*(a-b));

theta = acos(x0);

T = cos(theta.'*(0:1:N-1));
out = (T*cf).';


