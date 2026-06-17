function out = local_evans(lambda,s,p)

    fun = @(x,y)(A(x,lambda,s,p))*y;
    sol = ode23s(fun,[-p.h,0],[0;1],s.options);
    temp = sol.y(:,end);
    out = p.rho_p*(lambda^2*temp(2)+p.g*p.k^2*temp(1));



% hold on;
% plot(sol.x,sol.y(1,:),'-k','LineWidth',2);
% plot(sol.x,sol.y(2,:),'-r','LineWidth',2);
