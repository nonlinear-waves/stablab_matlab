function d = right_manifold(nu_left,nu_right,scale)

    t_start = tic; % start timer
    
    pie = iv('pi'); % rigorous pi
    
    M = 8;
    N = 2*M+1; % we set N = 2*M+1 so that we can use M in our induction proof on the bound on $\phi_{m,n}$.
    
    if N < 4
       error('Hypothesis of error bounds not met.'); 
    end 

    %
    % Constants
    %
    
    % reparameterize with alpha = sqrt(4*nu-1)
    temp = iv(maxi(inf(4*nu_left-1),0));
    alpha_left = sqrt(temp);
    alpha_right = sqrt(4*nu_right-1);
        
    % Determine a bound on the region where the coefficients of the series solution
    % are analytic in $\alpha$.
    imag_bound = 10^5;
    for m = 0:N
        for n = 0:N
            
            a = iv(-m^2 + 2*m*n - n^2 + 1);
            b = iv(-2*m^2 + 2*m + 2*n^2 - 2*n)*1j;
            c = iv(m^2 + 2*m*n - 2*m + n^2 - 2*n + 1);

            if a == 0
                imag_bound = min(imag_bound,inf(abs(imag(c/b))));
            else
                rt1 = (-b+sqrt(b*b-4*a*c))/(2*a);
                rt2 = (-b+sqrt(b*b-4*a*c))/(2*a);
                imag_bound = min(imag_bound,min(inf(abs(imag(rt1))),inf(abs(imag(rt2)))));
            end
       
        end
    end
    
    % Find $\rho$ such that the coefficients are analytic in alpha on and
    % insdie the ellipse $E_{\rho}$.
    alpha_diff = alpha_right-alpha_left;
    rho_max = 2*imag_bound/alpha_diff+sqrt(16*imag_bound^2/alpha_diff^2+4)/2;
    rho = iv(inf((2*rho_max)/3));
    
    % get the stadium for analytic interpolation error bound
    theta = linspace(0,sup(2*pie),1001);
    theta = iv(theta(1:end-1),theta(2:end));
    E_rho = (alpha_left+alpha_right)/2+(alpha_right-alpha_left)*(rho*exp(1i*theta)+exp(-1i*theta)/rho)/2;
    
    %
    % get a bound on |coefficients| on stadium
    %
    
    [phi_bound,w_bound,u_bound,z_bound] = coeff_right(E_rho,N,scale);
    
    M_rho = 1e-10;
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            if n > m
               continue 
            end
            t_phi = sup(abs(phi_bound(m+1,n+1)));
            t_w = sup(abs(w_bound(m+1,n+1)));
            t_u = sup(abs(u_bound(m+1,n+1)));
            t_z = sup(abs(z_bound(m+1,n+1)));
            M_rho=maxi([M_rho,t_phi,t_w,t_u,t_z]);
        end
    end
    
    
    %
    % find out how many Chebyshev nodes are needed for the desired error bound
    %
    
    D_rho = (rho+1/rho)/2-1;
    L_rho = pie*sqrt(rho^2+1/rho^2);

    M = 1; % The number of Chebyshev nodes
    eta = log(rho);
    interpolation_error = iv(1);
    while sup(interpolation_error) > 1e-17
        M = M+1;
        interpolation_error = M_rho*L_rho/(pie*D_rho*sinh(eta*(M+1)));
    end
    interpolation_error = iv(sup(interpolation_error));
    interp_err = iv(-interpolation_error,interpolation_error);
    
    %
    % interpolation nodes
    %
    
    theta = (2*iv(0:1:M-1)+1)*pie/(2*M);
    alpha_tilde = cos(theta);
    alpha_vals = (alpha_left+alpha_right)/2+(alpha_left-alpha_right)*alpha_tilde/2;
    
    %
    % Compute the coefficients on the interpolation nodes
    %
    
    [phi,w,u,z] = coeff_right(alpha_vals,N,scale);
    
    Tcf = get_Tcf(M);

    cf_phi = iv(zeros(size(phi)));
    cf_w = iv(zeros(size(w)));
    cf_u = iv(zeros(size(u)));
    cf_z = iv(zeros(size(z)));
    
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            cf_phi(m+1,n+1,:) = Tcf*squeeze(phi(m+1,n+1,:));
            cf_w(m+1,n+1,:) = Tcf*squeeze(w(m+1,n+1,:));
            cf_u(m+1,n+1,:) = Tcf*squeeze(u(m+1,n+1,:));
            cf_z(m+1,n+1,:) = Tcf*squeeze(z(m+1,n+1,:));
            
            cf_phi(m+1,n+1,1) = cf_phi(m+1,n+1,1) + interp_err;
            cf_w(m+1,n+1,1) = cf_w(m+1,n+1,1) + interp_err;
            cf_u(m+1,n+1,1) = cf_u(m+1,n+1,1) + interp_err;
            cf_z(m+1,n+1,1) = cf_z(m+1,n+1,1) + interp_err;
        end
    end
    
    
    % Prepare to evaluate the maximum of the coefficients over alpha
    % interval.
    theta_wide = iv(linspace(0,sup(2*pie),1001));
    theta_wide = iv(theta_wide(1:end-1),theta_wide(2:end));
    T_eval = cos(theta_wide.'*(0:1:M-1));
    
    % Determine C for $\phi$ in the inductive proof.
    C = 0;
    kappa = (1+alpha_right^2)/4;
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            if m^2+n^2 == 0
               continue 
            end
            max_phi = maxi(sup(abs(T_eval*squeeze(cf_phi(m+1,n+1,:)))));
            term1 = (kappa*max_phi)^(1/iv(m+n));
            C = maxi([C,sup(term1)]);
        end
    end
    
    C = iv(C);
    C0 = 4/(1+alpha_left^2);
    
    % Determine K, K0 for |w_{m,n}|\leq K0*K^(m+n) in the inductive proof.
    K = 1.05*C;
    a = 2/sqrt(1+alpha_left^2);
    
    % take K0 big enough that the induction proof holds
    K0 = sup(a*C0*N*(C/K)^N)+1e-6;
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            max_w = iv(maxi(sup(abs(T_eval*squeeze(cf_w(m+1,n+1,:))))));
            K0 = maxi([K0,max_w/K^(m+n)]);
        end
    end
    K0 = iv(K0);
    
    % Determine bounds on L, L0 such that |u_{m,n}| \leq L0*L^(m+n)
    L = 1.5*K;
    L_constant = K0*(1+alpha_right^2)^2/(iv(8)*N^2*(1-iv(K)/L)^2);
    if sup(L_constant) > 1
       error('Proof failed');
    end
    
    % Now take L_0 big enough to show the base cases of the induction proof
    L0 = iv(abs(u(1,1)));
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            max_u = iv(maxi(sup(abs(T_eval*squeeze(cf_u(m+1,n+1,:))))));
            L0 = maxi([L0,max_u/L^(m+n)]);
        end
    end
    L0 = iv(L0);
    
    % Determine R, R0 such that |z_{m,n}|\leq R0*R^(m+n) in the inductive proof.
    R = 1.05*L;
    R0 = sup(a*L0*N*(R/L)^N)+1e-6;
    for N_index = 0:N
        for m = 0:N_index
            n = N_index-m;
            max_z = iv(maxi(sup(abs(T_eval*squeeze(cf_z(m+1,n+1,:))))));
            R0 = maxi([R0,max_z/R^(m+n)]);
        end
    end
    R0 = iv(R0);
    
    %
    % combine the output into a structure
    %
    
    % Bound four terms with the same geometric series
    d.C = R;
    d.C0 = iv(maxi([sup(C0),sup(K0),sup(R0),sup(L0)]));
    
    d.N = N;
    
    % record the coefficients
    d.phi = phi;
    d.w = w;
    d.u = u;
    d.z = z;
    
    d.cf_phi = cf_phi;
    d.cf_w = cf_w;
    d.cf_u = cf_u;
    d.cf_z = cf_z;
    
    % record the compute time
    d.run_time = toc(t_start);
    
    
    
    