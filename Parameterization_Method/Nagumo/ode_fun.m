function out = ode_fun(x,y,P,n,lam)

% triple convolution
sm = 0;
for k = 0:n
    for l = 0:k
        if local_delta(n,k,l) == 1             
            temp1 = P{n-k+1}.fun(x);
            temp2 = P{k-l+1}.fun(x);
            temp3 = P{l+1}.fun(x);
            sm = sm + temp1(1)*temp2(1)*temp3(1);
        end
    end
end

B = [0, 1; 
    1+n*lam-6*sech(x)^2, 0];

out = B*y+[0;-sm];


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


