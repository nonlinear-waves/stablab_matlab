

This is a small Matlab script to calculate fronts for the following system:
> U' + cU - cz = 0  
> eps*V' + cV - cz = 0  
> z' = -(1/c)(1-V)*F(hU + (1-h)V)  

for eps=0. So, below, you should always use eps=0.

This system was produced by integration of
> U" + cU' + (1-V)F(W) = 0  
> eps*V" + cV' + (1-V)F(W) = 0  

where W = hU + (1-h)V.

This is just a very simple change of variable from the equation in the "older" paper

It shoots from near (1,1) and adjusts c to land at (0,0).

Usage:
-
    [x,u,v] = integrated_find_c(eps,h,Z,sigma,w_star);
Or with some fairly arbitrary numbers:

    [x,u,v] = integrated_find_c(0,0.3,6,0.25,1e-3);

Then you can do things like:

    plot(u,v);
    plot(x,[u v]);

