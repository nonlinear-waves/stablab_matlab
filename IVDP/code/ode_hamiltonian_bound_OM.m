function out = ode_hamiltonian_bound_OM(~,v,pars)

if norm(v(1:2))< 4

    x = v(1);
    y = v(2);
    p = v(3);
    q = v(4);

    d = fg(x,y,pars);

    out = [d.f+pars.sigma1^2*p;
        d.g+pars.sigma2^2*q;
        -d.f_x*p-d.g_x*q+4*pars.a*x*pars.eps;
        -d.f_y*p-d.g_y*q];

else
   
    out = zeros(4,1);
    
end











