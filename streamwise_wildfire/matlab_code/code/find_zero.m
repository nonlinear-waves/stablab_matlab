function c = find_zero(a,b,fun)

fa = fun(a);
fb = fun(b);

if fa*fb >= 0
    error('A zero may not be straddled.');
end


while abs(b-a) > 1e-12
    c = (a+b)/2;
    fc = fun(c);
    if fa*fc > 0
        a = c;
    else
        b = c;
    end
end





