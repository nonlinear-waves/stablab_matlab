function U_n = finite_diff_advance(U_n,U_o,K,H,p,tol,fd_F,fd_jac, ...
    bc_L_fun,bc_R_fun,bc_L_jac_fun,bc_R_jac_fun) 

old = U_n;

diff = tol+1;
while diff > tol

    % evaluate the system of finite difference equations
    Fout = fd_F(U_n,U_o,K,H,p);
    
    % evaluate the Jacobian of the system of finite difference equations
    Jout = fd_jac(U_n,U_o,K,H,p);
    
  
    % system size
    n = size(Fout,1);
    % number of interior nodes
    N = size(Fout,2);

    % Make the system of finite difference equations be in vector form
    Fval = reshape(Fout,n*N,1);

    
    % Get the Jacobian of Fval
    % This is slow and should eventually be changed for speed up!
    Jac = sparse(n*(N+2),n*(N+2));
    for j = 1:N % corresponds to node location
        for k = 1:3 % corresponds to j-1, j, j+1
            for l = 1:n % corresponds to system equation
                for r = 1:n % corresponds to system variable
%                     
%                     fprintf('\n(l,r,k) = (%4g,%4g,%4g)',l,r,k);
%                     Jout{l}{r}{k}
                    
                    if Jout{l}{r}{k} == 0
                        Jac(j*n+l,(j-1)*n+(k-1)*n+r) = 0; 
                    else                    
                        if length(Jout{l}{r}{k}) == 1
                            Jac(j*n+l,(j-1)*n+(k-1)*n+r) = Jout{l}{r}{k};
                        else
                            Jac(j*n+l,(j-1)*n+(k-1)*n+r) = Jout{l}{r}{k}(j);
                        end
                    end
                end
            end
        end
    end

    % left boundary conditions
    bc_L = bc_L_fun(U_n,U_o,K,H,p);
    % right boundary conditions
    bc_R = bc_R_fun(U_n,U_o,K,H,p);
    % Jacobian for left boundary conditions
    bc_L_jac = bc_L_jac_fun(U_n,U_o,K,H,p);
    % Jacobian for right boundary conditions
    bc_R_jac = bc_R_jac_fun(U_n,U_o,K,H,p);

    % system equations (adding boundary conditions)
    F = [bc_L;Fval;bc_R];
    % system Jacobian (adding Jacobian of boundary conditions)
    Jac(1:size(bc_L_jac,1),1:size(bc_L_jac,2)) = bc_L_jac;
    Jac(end-size(bc_R_jac,1)+1:end,end-size(bc_R_jac,2)+1:end) = bc_R_jac;
    
%     mx = max(max(abs(Jac)))

%     (full(Jac))
%     dt = det(Jac)
    
%     STOP
   
    % Newton solve
    U_n = reshape(reshape(U_n,n*(N+2),1)-Jac\F,n,N+2);
    % check solution convergence
    diff = max(norm(U_n-old));
    % update Newton iteration step
    old = U_n;

end

