function [un,u0,u1] = fd_Schrodinger(p,P,T,Del_x,L2,constant,theta1,theta2,lam1,stats)


Del_t = Del_x^2*constant;
xpoints = round(2*L2/Del_x);
tpoints = round(T/Del_t);

% Dependant params.
xgrid = linspace(-L2,L2,xpoints+1);
tgrid = linspace(0,T,tpoints+1);
K = tgrid(2) - tgrid(1);
H = xgrid(2) - xgrid(1);

% grid for the variable used in the homological equations 
z1 = theta1+1i*theta2;
z2 = conj(z1);
z3 = z1*exp(lam1*T);
z4 = conj(z3);
z_grid = z1*exp(lam1*tgrid);
z_grid2 = conj(z_grid);

% FD vars

tol = 1e-10;

%% form u0 and u1.

temp = evaluate_homological_equations_dim2(P,z1,z2,xgrid);
u0 = real(temp([1,3],:));
temp = evaluate_homological_equations_dim2(P,z3,z4,xgrid);
u1 = real(temp([1,3],:));

%% Time evolution finite difference progression.

un = real(u0);
bc_L_jac_fun = @(U_n,U_o,K,H,p)eye(2);
bc_R_jac_fun = @(U_n,U_o,K,H,p)eye(2);
for ind = 2:tpoints
    
    % Update the boundary values.
    z1 = z_grid(ind);
    z2 = z_grid2(ind);
        
    temp = evaluate_homological_equations_dim2(P,z1,z2,-L2);
    bc_L = real(temp([1,3],:));
    temp = evaluate_homological_equations_dim2(P,z1,z2,L2);
    bc_R = real(temp([1,3],:));

    % Create boundary functions
    bc_L_fun = @(U_n,U_o,K,H,p)(U_n(:,1)-bc_L);
    bc_R_fun = @(U_n,U_o,K,H,p)(U_n(:,end)-bc_R);
    
    % Update the initial guess.
    uo = real(un);
    un = finite_diff_advance(un,uo,K,H,p,tol,@fd_F,@fd_jac,bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun);
    
    if strcmp(stats,'on')
        clf;
        hold on;
        plot(xgrid, real(u0),'-r','LineWidth',2);
        plot(xgrid, real(u1),'-g','LineWidth',2);
        plot(xgrid, real(un),'-k','LineWidth',2);
        drawnow;
    end
    
end

if strcmp(stats,'on')
    fprintf('\n\n');
    fprintf(num2str(max(max(abs(real(un)-real(u1))))));
    fprintf('\n\n');
end
