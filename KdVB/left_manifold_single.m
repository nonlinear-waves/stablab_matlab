function d = left_manifold_single(nu,theta_0,r_newton)

    t_start = tic;

    %
    % Constants
    %
    
    theta_0_wide = theta_0+iv(-r_newton,r_newton);
    
    % Define some terms based on if the computation is rigorous.
    if isa(nu,'intval')
       ivfun = @(x)iv(x);
       rig = 1; % boolean that describes when the code is rigorous or not
       pie = iv('pi');
    else
        ivfun = @(x)x;
        rig = 0;
        pie = pi;
    end
    
    % check that |theta_0| < 1
    if sup(abs(theta_0)) >= 1
       error('By assumption, |\theta_0|<1.'); 
    end

    % Parameter mu 
    mu = (-1/(2*nu))*(1-sqrt(1+4*nu));
    
    q = 3; % Hypothesis: |a_n| \leq C_0C^n/(n+1)^q, (q >= 2 required, q = 3 efficient)
    C_0 = ivfun(75); % $C_0 \geq 1$ a choice. (Note that $0<|\mu|<1$.) 
    C = ivfun(1)/100; 
    
    % Take N large enough for the inductive hypothesis to hold
    %{
        
    %}
    N = 10; % N >= 2 required in assumptions about the truncation error
    xi_q = pie^2/6-1;
    while sup(2^q*C_0*xi_q*(ivfun(N+1)/(N+2))^q/(ivfun(N-1)*(N+1-N*mu))) > 1
        N = N+1;
    end

    % ---------------------------------------------------------------------
    % Now compute enough terms to verify the inductive hypothesis holds.
    % ---------------------------------------------------------------------
    
    % initialize phi
    phi = ivfun(zeros(1000,1));
    phi(1,:) = ivfun(1); % phi_0
    phi(2,:) = ivfun(-1); % phi_1
    
    % initialze w = phi'
    w = ivfun(zeros(1000,1));
    w(2,:) = -mu; % w_1, (w_0 = 0)

    % compute the rest of the terms 
    for n = 2:N

        % convolution $ \sum_{k=1}^{n-1} \phi_k \phi_{n-k} $
        sm = phi(2:n).'*flipud(phi(2:n));

        % phi_n = \frac{1}{2(n^2-1+n(1-n)\mu)} \sum_{k=1}^{n-1} \phi_k \phi_{n-k}$
        phi(n+1) = sm./(2*(ivfun(n)^2-1+ ivfun(n)*(1-n)*mu));
        w(n+1) = n*mu.*phi(n+1);
        
        % gurantee C is big enough for the inductive hypotheses to hold
        C = ivfun(maxi([sup(C); ...
            sup( (abs(phi(n+1))*(n+1)/C_0)^(1/iv(n)) ); ...
            sup( (abs(w(n+1))*(n+1)/C_0)^(1/iv(n)) )]));

    end
    
    % The trucation error requires that C < 1
    if sup(C) >= 1
       error('C is too big.'); 
    end
    
    % record C and C_0
    d.C = C;
    d.C_0 = C_0;
  
    % truncation error bound
    truncation_error = C_0*C^(N+1)/((1-C)*ivfun(N+2)^q);

    % compute more terms until the desired truncation error is met.
    while sup(truncation_error) > 1e-18
        
        n = n+1;
        % convolution $ \sum_{k=1}^{n-1} \phi_k \phi_{n-k} $
        sm = phi(2:n).'*flipud(phi(2:n));

        % phi_n = \frac{1}{2(n^2-1+n(1-n)\mu)} \sum_{k=1}^{n-1} \phi_k \phi_{n-k}$
        phi(n+1) = sm./(2*(ivfun(n)^2-1+ ivfun(n)*(1-n)*mu));
        w(n+1) = n*mu.*phi(n+1);
        
        truncation_error = C_0*C^(n+1)/((1-C)*ivfun(n+2)^q);
        
    end
    
    % record terms
    d.N = n;
    phi = phi(1:d.N+1);
    w = w(1:d.N+1);
    d.truncation_error = truncation_error;
    d.der_truncation_error = C_0*C^(N)/((1-C)*ivfun(N+2)^(q-1));
    d.der_der_truncation_error = C_0*C^(N-1)/((1-C)*ivfun(N+2)^(q-2));
    
    
    % compute phi, w, and theta derivatives at $x = 0$
    phi_0 = 0;
    w_0 = 0;
    phi_theta = 0;
    w_theta = 0;
    phi_theta_theta = 0;
    w_theta_theta = 0;
    for n = 0:length(phi)-1
       phi_0 = phi_0 + phi(n+1)*theta_0^n;
       w_0 = w_0 + w(n+1)*theta_0^n;
    end
    for n = 1:length(phi)-1
       phi_theta = phi_theta + phi(n+1)*n*theta_0^(n-1);
       w_theta = w_theta + w(n+1)*(n)*theta_0^(n-1);
    end
    for n = 2:length(phi)-2
       phi_theta_theta = phi_theta_theta + phi(n+1)*(n)*(n-1)*theta_0_wide^(n-2);
       w_theta_theta = w_theta_theta + w(n+1)*(n)*(n-1)*theta_0_wide^(n-2);  
    end

    % add on the truncation error if the computation is rigorous
    if rig == 1
        
        phi_0 = phi_0 + iv(-d.truncation_error,d.truncation_error);
        w_0 = w_0 + iv(-d.truncation_error,d.truncation_error);
        
        phi_theta = phi_theta + iv(-d.der_truncation_error,d.der_truncation_error);
        w_theta = w_theta + iv(-d.der_truncation_error,d.der_truncation_error);
        
        phi_theta_theta = phi_theta_theta + iv(-d.der_der_truncation_error,d.der_der_truncation_error);
        w_theta_theta = w_theta_theta + iv(-d.der_der_truncation_error,d.der_der_truncation_error);

    end

    
    % record quantities
    
    d.phi = phi;
    d.w = w;
    
    d.phi_0 = phi_0;
    d.w_0 = w_0;
    d.phi_theta = phi_theta;
    d.w_theta = w_theta;
    d.phi_theta_theta = phi_theta_theta;
    d.w_theta_theta = w_theta_theta;
   
    % note the run time in seconds
    d.run_time = toc(t_start);

end























