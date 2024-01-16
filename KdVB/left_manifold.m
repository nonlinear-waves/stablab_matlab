function d = left_manifold(nu_left,nu_right,bar_theta_0)

    t_start = tic;
    
    %
    % Constants
    %

    pie = iv('pi'); % rigorous enclosure of $\pi$
    
    % 
    % Parameter bounds
    %

    mu_right = (-1/(2*nu_left))*(1-sqrt(1+4*nu_left));  % right bound on mu
    mu_right = iv(sup(mu_right));
    mu_left = (-1/(2*nu_right))*(1-sqrt(1+4*nu_right)); % left bound on mu
    mu_left = iv(inf(mu_left));
    
    % mu_tilde \in [-1,1], but we extend the domain for analytic
    % interpolation
    mu_tilde_left_bound = -(mu_left+mu_right)/(mu_right-mu_left); % guarantee mu > 0
    mu_tilde_right_bound = (2-(mu_left+mu_right))/(mu_right-mu_left); % gurantee mu < 1
    
    % guarantee that 0 < mu < 1
    max_mu_tilde = min(inf(abs(mu_tilde_left_bound+1)),inf(abs(mu_tilde_right_bound)));
    
    %
    % radius of the stadium for Chebyshev interpolation
    %

    % make rho close to optimal
    rho = iv(1);
    if sup((rho+1/rho)/2) >= max_mu_tilde
       error('This step faild.') 
    end
    
    while sup((rho+1/rho)/2)< max_mu_tilde
       rho_old = rho;
       rho = iv(sup(1.05*rho_old)); 
    end
    rho = rho_old;
    
    %
    % Chebyshev interpolation error
    %

    % find out how many Chebyshev nodes are needed for the desired error bound
    D_rho = (rho+1/rho)/2-1;
    L_rho = pie*sqrt(rho^2+1/rho^2);
    
    %
    % Finite series
    %

    N = ceil((log(0.5)+log(1e-17))/log(0.5)); % number of terms in series to keep

 
    if sup(bar_theta_0) > 1
        error('bar_theta_0 <= 1 is violated');
    end
    if N <= 2
       error('N > 2 is violated'); 
    end
    
    M = 1; % The number of Chebyshev nodes
    eta = log(rho);
    interpolation_error = iv(1);
    ratio = iv(sup(abs(bar_theta_0/2)));
    
    M_rho = 1+(iv(N)*iv(N-1)*iv(N-2)/bar_theta_0^2)*((1-ratio^N)/(1-ratio));
    while sup(interpolation_error) > 1e-16
        M = M+1;
        interpolation_error = M_rho*L_rho/(pie*D_rho*sinh(eta*(M+1)));
    end
    
    %
    % interpolation nodes
    %
    
    theta = (2*iv(0:1:M-1)+1)*pie/(2*M);
    mu_tilde = cos(theta);
    mu = (mu_left+mu_right)/2+(mu_right-mu_left)*mu_tilde/2;

    % initialize phi
    phi = iv(zeros(N+1,length(mu)));
    phi(1,:) = iv(1);
    phi(2,:) = iv(-1);

    for n = 2:N

        % convolution $ \sum_{k=1}^{n-1} \phi_k \phi_{n-k} $
        sm = sum(phi(2:n,:).*flipud(phi(2:n,:)),1);

        % phi_n = \frac{1}{2(n^2-1+n(1-n)\mu)} \sum_{k=1}^{n-1} \phi_k \phi_{n-k}$
        phi(n+1,:) = sm./(2*(iv(n)^2-1+ iv(n)*(1-n)*mu));

    end
    
    % compute the finite series involving phi
    phi_sum = iv(zeros(M,1));
    phi_theta_sum = iv(zeros(M,1));
    phi_theta_theta_sum = iv(zeros(M,1));
    phi_theta_theta_theta_sum = iv(zeros(M,1));
    for n = fliplr(0:N)
        phi_sum = phi_sum + bar_theta_0^n*phi(n+1,:).';
    end
    
    for n = fliplr(1:N)
        phi_theta_sum = phi_theta_sum + n*bar_theta_0^(n-1)*phi(n+1,:).';
    end
    
    for n = fliplr(2:N)
        phi_theta_theta_sum = phi_theta_theta_sum + n*(n-1)*bar_theta_0^(n-2)*phi(n+1,:).';
    end
     
    for n = fliplr(3:N)
        phi_theta_theta_theta_sum = phi_theta_theta_theta_sum + ...
            n*(n-1)*(n-2)*bar_theta_0^(n-3)*phi(n+1,:).';
    end
    
    % get the transformation matrix from function values to Chebyshev
    % coefficients
    Tcf = get_Tcf(M);
    
    % compute the Chebyshev coefficients for phi
    cf_phi = Tcf*phi_sum;
    cf_phi_theta = Tcf*phi_theta_sum;
    cf_phi_theta_theta = Tcf*phi_theta_theta_sum;
    cf_phi_theta_theta_theta = Tcf*phi_theta_theta_theta_sum;
    
    % add on the enclosure of the interpolation error for phi
    err = iv(-interpolation_error,interpolation_error);
    cf_phi(1) = cf_phi(1)+ err;
    cf_phi_theta(1) = cf_phi_theta(1) + err;
    cf_phi_theta_theta(1) = cf_phi_theta_theta(1) + err;
    cf_phi_theta_theta_theta(1) = cf_phi_theta_theta_theta(1) + err;
    
    % add on the enclosure of the truncation error for phi and theta
    % derivatives
    temp = f_bound(bar_theta_0/2,N);
    cf_phi(1) = cf_phi(1)+ + iv(-temp(1),temp(1));
    cf_phi_theta(1) = cf_phi_theta(1) + iv(-temp(2),temp(2));
    cf_phi_theta_theta(1) = cf_phi_theta_theta(1) + iv(-temp(3),temp(3));
    cf_phi_theta_theta_theta(1) = cf_phi_theta_theta_theta(1)+ iv(-temp(4),temp(4));
    
    % get Chebyshev coefficients for mu
    cf_mu = Tcf*mu.';
    
    % get Chebyshev coefficients for w and its theta derivatives
    cf_w = poly_mult(cf_mu,bar_theta_0*cf_phi_theta);

    cf_w_theta = poly_add(poly_mult(cf_mu,cf_phi_theta), ...
        poly_mult(cf_mu,bar_theta_0*cf_phi_theta_theta));
    
    cf_w_theta_theta = poly_add(poly_mult(2*cf_mu,cf_phi_theta_theta), ...
        poly_mult(cf_mu,bar_theta_0*cf_phi_theta_theta_theta));
 
    % record the varios quantities
    d.cf_phi = cf_phi;
    d.cf_phi_theta = cf_phi_theta;
    d.cf_phi_theta_theta = cf_phi_theta_theta;
    d.cf_phi_theta_theta_theta = cf_phi_theta_theta_theta;
    
    d.cf_w = cf_w;
    d.cf_w_theta = cf_w_theta;
    d.cf_w_theta_theta = cf_w_theta_theta;
    
    d.nu_left = nu_left;
    d.nu_right = nu_right;
    d.mu_left = mu_left;
    d.mu_right = mu_right;
    d.mu_tilde_left_bound = mu_tilde_left_bound;
    d.mu_tilde_right_bound = mu_tilde_right_bound;
    d.rho = rho;
    d.N = N;
    d.M = M;
    d.interpolation_error = interpolation_error;
    d.theta = theta;
    d.mu_tilde = mu_tilde;
    d.mu = mu;
    d.nu = (1-mu)./mu.^2;
    d.phi = phi;
    d.bar_theta_0 = bar_theta_0;
    d.run_time = toc(t_start);
    
    
end

function out = f_bound(x,N)

    % Contains $f(x):= x^(N+1)/(1-x)$, and its first three derivatives
    out = [x^(N+1)/(1-x); 
            (N+1)*x^N/(1-x)+x^(N+1)/(1-x)^2;
            N*(N+1)*x^(N-1)/(1-x)+2*(N+1)*x^N/(1-x)^2+2*x^(N+1)/(1-x)^3;
            N*(N+1)*(N-1)*x^(N-2)/(1-x)+3*N*(N+1)*x^(N-1)/(1-x)^2 + ...
            6*(N+1)*x^N/(1-x)^3+6*x^(N+1)/(1-x)^4];
        
end
 
    
    
    
    
    




