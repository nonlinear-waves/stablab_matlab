function [U_n,u0,u1] = fd_Gray_Scott(p,P,T,H,constant,sigma0,L_fd,lam,stats)


%% dependent constants and values
K = H^2*constant;
sigma1 = sigma0*exp(lam*T);

tpoints = round(T/K);
xpoints = round(2*L_fd/H);
xgrid = linspace(-L_fd,L_fd, xpoints);
tgrid = linspace(0,T,tpoints);
% Set up grid
K = tgrid(2)-tgrid(1);
H = xgrid(2)-xgrid(1);
lambdagrid = sigma0*exp(lam*tgrid);

%% Produce sigma values.

if strcmp(stats,'on')
    figure;
    hold on;
end

u0 = evaluate_homological_equations(P,sigma0,xgrid);
u1 = evaluate_homological_equations(P,sigma1,xgrid);

% plot(xgrid,u0(1,:),'-k','LineWidth',2);
% plot(xgrid,u1(1,:),'-g','LineWidth',2);
% drawnow;


%% Time evolution

U_n = [u0(1,:);u0(3,:)];
U_o = U_n;
tol = 1e-8;

bc_L_jac_fun = @(U_n,U_o,K,H,p) [1 0; 0 1];
bc_R_jac_fun = @(U_n,U_o,K,H,p) [1 0; 0 1];
for j = 1:length(tgrid)

    % Update the boundary values.
    sigma_n = lambdagrid(j);
    
    bc_L = evaluate_homological_equations(P,sigma_n,xgrid(1));
    bc_R = evaluate_homological_equations(P,sigma_n,xgrid(end));
   
    % Create boundary functions
    bc_L_fun = @(U_n,U_o,K,H,p)[U_n(1,1)-bc_L(1);U_n(2,2)-bc_L(3)];
    bc_R_fun = @(U_n,U_o,K,H,p)[U_n(end-1,end-1)-bc_R(1);U_n(end,end)-bc_R(3)];
    
    % Update finite difference
    U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,@fd_F,@fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);

    % Plot current 
    if strcmp(stats,'on')
%         j
        clf;
        hold on;
        plot(xgrid,u1([1,3],:),'-g','LineWidth',2);
        plot(xgrid,u0([1,3],:),'-k','LineWidth',2);
        plot(xgrid,U_n,'--r','LineWidth',2);
        drawnow;
    end
    
    U_o = U_n;

end