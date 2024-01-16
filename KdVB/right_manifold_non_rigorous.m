function [fun,d,fun_der] = right_manifold_non_rigorous(nu,scale)

    %
    % Constants
    %

    % transform the parameter nu to the parameter alpha
    alpha = sqrt(4*nu-1);
    
    % Take $N$ large enough to guarantee inductive hypothesis holds.
    max_N = 30;
    
    % Initialize coefficients of series
    phi = zeros(max_N+1,max_N+1);
    w = phi;
    u = phi;
    z = phi;
    
    % fixed point
    phi(1,1) = -1;
    w(1,1) = 0;
    u(1,1) = 1;

    % eigenvectors of the profile Jacobian at positive infinity
    phi(2,1) = -scale;
    phi(1,2) = -scale;

    mu1 = -1/(2*nu)-1i*sqrt(4*nu-1)./(2*nu);
    mu2 = conj(mu1);

    w(2,1) = -scale*mu1;
    w(1,2) = -scale*mu2;
    
    u(2,1) = w(2,1)/(2*mu1^2);
    u(1,2) = w(1,2)/(2*mu2^2);
    z(2,1) = mu1*u(2,1);
    z(1,2) = mu2*u(1,2);
    
    % Now obtain the coefficients of the series and obtain C so the lemma holds.
    C1 = (1+alpha^2)/2;
    C2 = 1+1i*alpha;
    C3 = 1-1i*alpha;
    C4 = 1+alpha^2;
    for N_ind = 2:max_N
       for n = 0:N_ind
            m = N_ind-n;

            % $\sum_{j = 0}^m \sum_{k = 0}^n \delta_{j,k} \phi_{j,k} \phi_{m-j,n-k}$
            sm = sum(diag(phi(1:m+1,1:n+1)*rot90(phi(1:m+1,1:n+1),2).'))-2*phi(1,1)*phi(m+1,n+1);

            % $\phi_{m,n}$ and $w_{m,n}$     
            Omega = C1/((m*C2+n*C3)^2-2*m*C2-2*n*C3+C4);
            phi(m+1,n+1) = Omega.*sm;
            w(m+1,n+1) = (m*mu1+n*mu2)*phi(m+1,n+1);

            % u and z
            sm_u = sum(diag(w(1:m+1,1:n+1)*rot90(u(1:m+1,1:n+1),2).'))-w(1,1)*u(m+1,n+1)-w(m+1,n+1)*u(1,1);    
            Omega_u = 0.5./((m*mu1+n*mu2).^2);
            u(m+1,n+1) = Omega_u.*sm_u;
            z(m+1,n+1) = (m*mu1+n*mu2)*u(m+1,n+1);
                
       end
    end
    
    fun = @(theta,r)real(eval(theta,r,phi,w,u,z));
    
    d.phi = phi;
    d.w = w;
    d.u = u;
    d.z = z;
    
    if nargout > 2
        fun_der = @(theta,r)real(eval_der(theta,r,phi,w,u,z));

    end

end
    
% evaluate the parameterization
function out = eval(theta,r,phi,w,u,z)

    theta_1 = r*(cos(theta)+1i*sin(theta));
    theta_2 = conj(theta_1);

    T1 = theta_1.^(0:1:size(phi,1)-1);
    T2 = theta_2.^(0:1:size(phi,2)-1);

    phi_val = T1*phi*T2.';
    w_val = T1*w*T2.';
    u_val = T1*u*T2.';
    z_val = T1*z*T2.';

    out = [phi_val;w_val;u_val;z_val];

end

% evaluate the derivative of the evaluation
function out = eval_der(theta,r,phi,w,u,z)

    theta_1 = r*(cos(theta)+1i*sin(theta));
    theta_2 = conj(theta_1);

    theta_1_der = r*(-sin(theta)+1i*cos(theta));
    theta_2_der = conj(theta_1_der);
    
    T1 = theta_1.^(0:1:size(phi,1)-1);
    T2 = theta_2.^(0:1:size(phi,2)-1);
    
    T1_der = (0:1:size(phi,1)-1).*theta_1.^[0,0:1:size(phi,1)-2].*theta_1_der;
    T2_der = (0:1:size(phi,2)-1).*theta_2.^[0,0:1:size(phi,2)-2].*theta_2_der;

  
    phi_val = T1_der*phi*T2.'+T1*phi*T2_der.';
    w_val = T1_der*w*T2.'+T1*w*T2_der.';
    u_val = T1_der*u*T2.'+T1*u*T2_der.';
    z_val = T1_der*z*T2.'+T1*z*T2_der.';

    out = [phi_val;w_val;u_val;z_val];

end









