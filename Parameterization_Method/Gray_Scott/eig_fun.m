function out = eig_fun(dom,sol)

out = zeros(length(dom),4);

for j = 1:length(dom)
    
    x = dom(j);

    if x < 0
        temp = deval(sol,x);
        out(j,:) = temp(1:4);
    else 
        temp = deval(sol,sol.x(1)+x);
        out(j,:) = temp(5:8);
    end


end




