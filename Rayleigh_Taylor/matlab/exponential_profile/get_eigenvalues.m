function lambda = get_eigenvalues(index,k,a,h)

eps = 1e-10;
gravitational_constant = 9.8;


x1 = -pi/h+2*pi*index/h+eps;
x2 = pi/h+2*pi*index/h-eps;

g1 = (4*k^2-a^2+x1^2)*tan(x1*h/2)+2*a*x1;
g2 = (4*k^2-a^2+x2^2)*tan(x2*h/2)+2*a*x2;

if g1*g2 >= 0
    error('root is not bracketed');
end

% bisection method
while abs((x2-x1)/x1) > 1e-12
    xmid = (x1+x2)/2;
    gmid = (4*k^2-a^2+xmid^2)*tan(xmid*h/2)+2*a*xmid;
    if sign(gmid) == sign(g2)
        x2 = xmid;
    else
        x1 = xmid;
    end
end

x = (x1+x2)/2;

% eigenvalue 
lambda = sqrt(4*k^2*a*gravitational_constant/(a^2+4*k^2+x^2));


