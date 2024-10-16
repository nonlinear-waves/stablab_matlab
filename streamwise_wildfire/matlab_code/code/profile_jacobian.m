function out = profile_jacobian(U,c,x,p)

[fv,df] = f(U(1),p);

out = [0,1; 
    -df, -(c-Wz(x,p))];



