function u0 = evaluate_homological_equations_dim2(P,z1,z2,xgrid)


u0 = P{1}{1}.fun(xgrid);
for i= 2:length(P)
   for j=1:i
       k = i-j+1;
            u0 = u0 + z1^(j-1)*z2^(k-1)*P{j}{k}.fun(xgrid);
    end
end

    

