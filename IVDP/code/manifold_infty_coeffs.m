function P = manifold_infty_coeffs(p,N,scale)

% Find the eigenvalues with positive real part and their eigenvectors
DF = hamiltonian_jacobian(zeros(4,1),p);
[V,D] = eigs(DF)
ind = find(diag(D)>0);
xi_1 = V(:,ind(1));
xi_2 = V(:,ind(2));

% make it so the results are consistent regardless of the size of
% epsilon_1,2
xi_1 = xi_1*((0.046675870782923 - 0.016168995936709i)/xi_1(1));
xi_2 = xi_2*((0.046675870782923 + 0.016168995936709i)/xi_2(1));

lambda_1 = D(ind(1),ind(1));
lambda_2 = D(ind(2),ind(2));

% initialize the coefficients
P = zeros(N+1, N+1, 4);
 
% The fixed point is zero, initialize the eigenvectors.
P(2, 1, :) = scale*xi_1;
P(1, 2, :) = scale*xi_2;

Id = eye(4);       
C = 0;
C0 = 10;
for  order = 2:N

     for m = 0:order
         
        n = order - m;

        A = DF - (n*lambda_1 + m*lambda_2)*Id;

        Sum1 = cauchy_product_3(P(:, :, 1), P(:, :, 1), P(:, :, 2), n, m);
        Sum2 = cauchy_product_3(P(:, :, 1), P(:, :, 2), P(:, :, 4), n, m);
        Sum3 = cauchy_product_3(P(:, :, 1), P(:, :, 1), P(:, :, 4), n, m);

       Rmn = [0;
              -2*p.a*Sum1;
              4*p.a*Sum2;
              2*p.a*Sum3];

        P(n+1, m+1, :) = A\Rmn;
        
        temp = abs((squeeze(P(n+1,m+1,:)))./C0).^(1/order);
        C = max([C;temp]);
     end

end

if C > 1
   error('Radius of convergence too big for how we use it.'); 
end

if C0*C^N > 1e-16
   error('N not big enough'); 
end
