function out = ode_fun(x,y,n,m,lam1,lam2,P,p)


% not designed for efficiency

ind = 1;

temp = P{1}{1}.fun(x);
u0 = temp(1);
v0 = temp(3);

u = zeros(n,m);
v = zeros(n,m);
sm = n+m;
for diag_number = 0:sm-1
    for j = 0:diag_number
        k = diag_number - j;
        temp = P{j+ind}{k+ind}.fun(x);
        u(j+ind,k+ind) = temp(1);
        v(j+ind,k+ind) = temp(3);
    end
end

M = [0,1,0,0;
    1-3*u0^2-v0^2, 0, (lam1*n+lam2*m)+2*p.nu-2*u0*v0, 0;
    0, 0, 0, 1;
    -(lam1*n+lam2*m)-2*u0*v0, 0, p.mu-3*v0^2-u0^2, 0];
    
N1 = 0;
N2 = 0;
for k = 0:n
    for j = 0:k
        for l = 0:m
            for q = 0:l

                if ~( (k^2+l^2)*((j-n)^2+(q-m)^2)*((k-n)^2+(l-m)^2+j^2+q^2)  == 0)

                    N1 = N1 + v(n-k+ind,m-l+ind)*v(k-j+ind,l-q+ind)*v(j+ind,q+ind)...
                        +v(n-k+ind,m-l+ind)*u(k-j+ind,l-q+ind)*u(j+ind,q+ind);
                    
                    N2 = N2 + u(n-k+ind,m-l+ind)*u(k-j+ind,l-q+ind)*u(j+ind,q+ind)...
                        +u(n-k+ind,m-l+ind)*v(k-j+ind,l-q+ind)*v(j+ind,q+ind);

                end
            end
        end
    end
end


out = M*y-[0; N2; 0; N1];
    
    