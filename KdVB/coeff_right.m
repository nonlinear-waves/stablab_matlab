function [phi,w,u,z] = coeff_right(alpha,max_N,scale)
%{
    [phi,w,u,z] = coeff_right(alpha,max_N,scale)

    Returns the coefficients of the series solution for the manifold at
    positive infinity. The input is alpha where nu = (alpha^2+1)/4, max_N
    which is the number of terms in the series we keep, and scale which is
    some number just barely smaller than 0.1. Returns the coefficients
    evaluated at alpha.
%}

% Rigorous representation of 1/2
one = iv(1);
half = one/2;

J = length(alpha);

% Initialize coefficients of series
phi = iv(zeros(max_N+1,max_N+1,J));
w = iv(zeros(max_N+1,max_N+1,J));
u = iv(zeros(max_N+1,max_N+1,J));
z = iv(zeros(max_N+1,max_N+1,J));

% fixed point
phi(1,1,:) = -one;
w(1,1,:) = iv(0);
u(1,1,:) = one;
z(1,1,:) = iv(0);

% eigenvalues of profile Jacobian at positive infinity
mu1 = (-2*(1-1i*alpha)./(1+alpha.^2)).';
mu2 = (-2*(1+1i*alpha)./(1+alpha.^2)).';


% nu in terms of alpha
nu = (alpha.^2+1).'/4;

% eigenvectors of the profile Jacobian at positive infinity
phi(2,1,:) = -scale;
phi(1,2,:) = -scale;    

w(2,1,:) = -mu1*scale;
w(1,2,:) = -mu2*scale;

u(2,1,:) = squeeze(w(2,1,:))./(2*mu1.^2);
u(1,2,:) = squeeze(w(1,2,:))./(2*mu2.^2);

z(2,1,:) = mu1.*squeeze(u(2,1,:));
z(1,2,:) = mu2.*squeeze(u(1,2,:));


% all the other coefficients
for N_ind = 2:max_N
    
    for n = 0:N_ind
        m = N_ind-n;
        
        % $\sum_{j = 0}^m \sum_{k = 0}^n \delta_{j,k} \phi_{j,k} \phi_{m-j,n-k}$
        sm = iv(zeros(1,1,J));
        for j = 0:m
            for k = 0:n
                if (j==0)&&(k==0)
                    continue
                end
                if (j==m)&&(k==n)
                    continue
                end
                sm = sm + phi(j+1,k+1,:).*phi(m-j+1,n-k+1,:);
            end
        end
        sm = squeeze(sm);
        
        % $\phi_{m,n}$ and $w_{m,n}$
        Omega = half./(nu.*(m*mu1+n*mu2).^2+m*mu1+n*mu2+1);
        
        
        
        phi(m+1,n+1,:) = Omega.*sm;
        
        w(m+1,n+1,:) = (m*mu1+n*mu2).*squeeze(phi(m+1,n+1,:));
        
        % $u_{m,n}$ and $z_{m,n}$

        sm_u = iv(zeros(1,1,J));
        for j = 0:m
            for k = 0:n
                if (j==0)&&(k==0)
                    continue
                end
                if (j==m)&&(k==n)
                    continue
                end
                sm_u = sm_u + w(j+1,k+1,:).*u(m-j+1,n-k+1,:);
            end
        end
       
        Omega_u = half./((m*mu1+n*mu2).^2);
        u(m+1,n+1,:) = Omega_u.*squeeze(sm_u);
        z(m+1,n+1,:) = (m*mu1+n*mu2).*squeeze(u(m+1,n+1,:));
        
    end
end






