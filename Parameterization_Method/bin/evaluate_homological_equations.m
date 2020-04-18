function u0 = evaluate_homological_equations(P,theta,xgrid)



% Plot the solution at theta
u0 = zeros(length(P{1}.fun(xgrid(1))),length(xgrid));
for j = 0:length(P)-1
    u0 = u0+theta^j*P{j+1}.fun(xgrid);
end
    

