function out = ode_fun(x,y,P,n,lam,p)

% triple convolution
sm = 0;
for k = 0:n
    for l = 0:k
        if local_delta(n,k,l) == 1             
            temp1 = P{n-k+1}.fun(x);
            temp2 = P{k-l+1}.fun(x);
            temp3 = P{l+1}.fun(x);
            sm = sm + temp1(1)*temp2(3)*temp3(3);            
        end
    end
end

temp = P{1}.fun(x);
a0 = temp(1);
b0 = temp(3);

B = [0, 1, 0, 0; 
    (lam*n+p.alpha+b0^2), 0, 2*b0*a0, 0;
    0, 0, 0, 1;
    -b0^2/p.gamma, 0, (lam*n+(1-2*b0*a0)/p.gamma), 0];

out = B*y+[0;sm;0;-sm/p.gamma];


% delta function
function out = local_delta(n,k,l)


if (k == 0) && (l == 0)
    out = 0;
elseif (k == n) && (l == n)
    out = 0;
elseif (k == n) && (l == 0)
    out = 0;
else
   out = 1;
end


