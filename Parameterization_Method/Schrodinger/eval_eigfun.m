function out = eval_eigfun(sol,x,L)


out = zeros(4,length(x));
for j = 1:length(x)
    if x(j) < 0
        temp = deval(sol,x(j));
        out(:,j) = temp(1:4);
    else
        temp = deval(sol,x(j)-L);
        out(:,j) = temp(5:8);
    end
end