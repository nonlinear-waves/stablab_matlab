function out = eval_rm(theta,r,dr)

% for convenience
N = dr.N;
C = dr.C;
C0 = dr.C0;

% Get terms that are needed to compute theta_{\pm}
rc = poly_mult(r,cos_of_poly(theta));
rs = poly_mult(r,sin_of_poly(theta));

rcon1 = poly_add(-rs,1i*rc);
rcon2 = poly_add(-rs,-1i*rc);
rcon3 = poly_add(-rc,-1i*rs);
rcon4 = poly_add(-rc,1i*rs);

%{

r_R_alpha = 0.109
theta1 = -0.0590 + 0.0918*1i
R = 0.0474
nu = 0.26


nu = 0.25

r_R_alpha = 0.35
theta_R_alpha = 4.569
rc = -0.0140
rs = -0.351
theta_1 = -0.0139-0.351i
 C = 0.414
 C_0 = 77.9
 N = 17
 R = 0.751
 err = 40.6
%}

% Get R for the geometric sum remainder bound
theta1 = poly_add(rc,1i*rs);
theta2 = poly_add(rc,-1i*rs);
R = iv(sup(C*sum(abs(theta1))));
if sup(R) >= 1
   C
   sum(abs(theta1))
   disp(R);
   error('R too big'); 
end

% compute truncation of series error
err = C0*(N+2)*R^(N+1)/(1-R)+C0*R^(N+2)/(1-R)^2;


% pre-compute powers of theta1 and theta2
theta1_m{N+1} = iv(0);
theta2_n{N+1} = iv(0);
theta1_m{1} = iv(1);
theta2_n{1} = iv(1);
for m = 1:N
    theta1_m{m+1} = poly_mult(theta1_m{m},theta1);
    theta2_n{m+1} = poly_mult(theta2_n{m},theta2);
end

% get evaluation of phi, w, u, z
y_phi = iv(zeros(2,1));
y_w = iv(zeros(2,1));
y_u = iv(zeros(2,1));
y_z = iv(zeros(2,1));

for N_index = fliplr(0:N)
   for m = 0:N_index
      n = N_index-m;
      
      temp = poly_mult(theta1_m{n+1},theta2_n{m+1});

      y_phi = poly_add(y_phi,poly_mult(squeeze(dr.cf_phi(m+1,n+1,:)),temp));
      y_w = poly_add(y_w,poly_mult(squeeze(dr.cf_w(m+1,n+1,:)),temp));
      y_u = poly_add(y_u,poly_mult(squeeze(dr.cf_u(m+1,n+1,:)),temp));
      y_z = poly_add(y_z,poly_mult(squeeze(dr.cf_z(m+1,n+1,:)),temp));
      
   end
end

% add on error
y_phi = real(y_phi);
y_phi(1) = y_phi(1)+iv(-err,err);
y_w = real(y_w);
y_w(1) = y_w(1)+iv(-err,err);
y_u = real(y_u);
y_u(1) = y_u(1)+iv(-err,err);
y_z = real(y_z);
y_z(1) = y_z(1)+iv(-err,err);

out.y_phi = clip_tail(y_phi);
out.y_w = clip_tail(y_w);
out.y_u = clip_tail(y_u);
out.y_z = clip_tail(y_z);

% get evaluation of phi, w, u, z first derivatives with respect to theta

y_phi_theta = iv(zeros(2,1));
y_w_theta = iv(zeros(2,1));
y_u_theta = iv(zeros(2,1));
y_z_theta = iv(zeros(2,1));


for M = 1:N
   for m = 0:M
      n = M-m;
      
      if n == 0
        temp1 = iv(zeros(2,1));  
      else
        temp1 = n*poly_mult(rcon1,poly_mult(theta1_m{n},theta2_n{m+1}));
      end
      
      if m == 0
        temp2 = iv(zeros(2,1));
      else
        temp2 = m*poly_mult(rcon2,poly_mult(theta1_m{n+1},theta2_n{m}));
      end
      temp = poly_add(temp1,temp2);
      
      y_phi_theta = poly_add(y_phi_theta,poly_mult(squeeze(dr.cf_phi(m+1,n+1,:)),temp));
      y_w_theta = poly_add(y_w_theta,poly_mult(squeeze(dr.cf_w(m+1,n+1,:)),temp));
      y_u_theta = poly_add(y_u_theta,poly_mult(squeeze(dr.cf_u(m+1,n+1,:)),temp));
      y_z_theta = poly_add(y_z_theta,poly_mult(squeeze(dr.cf_z(m+1,n+1,:)),temp));

   end
end

% add on error
y_phi_theta = real(y_phi_theta);
y_phi_theta(1) = y_phi_theta(1)+iv(-err,err);
y_w_theta = real(y_w_theta);
y_w_theta(1) = y_w_theta(1)+iv(-err,err);
y_u_theta = real(y_u_theta);
y_u_theta(1) = y_u_theta(1)+iv(-err,err);
y_z_theta = real(y_z_theta);
y_z_theta(1) = y_z_theta(1)+iv(-err,err);

out.y_phi_theta = clip_tail(y_phi_theta);
out.y_w_theta = clip_tail(y_w_theta);
out.y_u_theta = clip_tail(y_u_theta);
out.y_z_theta = clip_tail(y_z_theta);

% get evaluation of phi, w, u, z second derivatives with respect to theta

y_phi_theta_theta = iv(zeros(2,1));
y_w_theta_theta = iv(zeros(2,1));
y_u_theta_theta = iv(zeros(2,1));
y_z_theta_theta = iv(zeros(2,1));

for M = 2:N
   for m = 0:M
      n = M-m;
      
      if (n==0)||(n==1)
        temp1 = iv(zeros(2,1));
      else
        temp1 = n*(n-1)*poly_mult(poly_mult(theta1_m{n-1}, ...
            theta2_n{m+1}),rcon1);
      end
      
      if (m==0)||(n==0)
        temp2 = iv(zeros(2,1));
      else
        temp2 =  -2*m*n*poly_mult(poly_mult(theta1_m{n}, ...
            theta2_n{m}),rs);
      end
      
      if (m==0)||(m==1)
        temp3 = iv(zeros(2,1));
      else
        temp3 = m*(m-1)*poly_mult(poly_mult(theta1_m{n+1}, ...
            theta2_n{m-1}),rcon2);
      end
      
      
      if (n==0)
        temp0 = iv(zeros(2,1));
      else
        temp0 = n*poly_mult(poly_mult(theta1_m{n},theta2_n{m+1}),rcon3);
      end
      
      if (m==0)
        temp4 = iv(zeros(2,1));
      else
        temp4 = m*poly_mult(poly_mult(theta1_m{n+1},theta2_n{m}),rcon4);
      end
      
    
      temp5 = poly_add(temp1,temp2);
      
      temp6 = poly_add(temp3,temp4);
      
      temp = poly_add(poly_add(temp0,temp5),temp6);
      
      y_phi_theta_theta = poly_add(y_phi_theta_theta, ...
          poly_mult(squeeze(dr.cf_phi(m+1,n+1,:)),temp));
      
        y_w_theta_theta = poly_add(y_w_theta_theta, ...
      poly_mult(squeeze(dr.cf_w(m+1,n+1,:)),temp));
      
        y_u_theta_theta = poly_add(y_u_theta_theta, ...
      poly_mult(squeeze(dr.cf_u(m+1,n+1,:)),temp));
      
        y_z_theta_theta = poly_add(y_z_theta_theta, ...
      poly_mult(squeeze(dr.cf_z(m+1,n+1,:)),temp));

   end
end

% add on error
y_phi_theta_theta = real(y_phi_theta_theta);
y_phi_theta_theta(1) = y_phi_theta_theta(1) + iv(-err,err);
y_w_theta_theta = real(y_w_theta_theta);
y_w_theta_theta(1) = y_w_theta_theta(1) + iv(-err,err);
y_u_theta_theta = real(y_u_theta_theta);
y_u_theta_theta(1) = y_u_theta_theta(1) + iv(-err,err);
y_z_theta_theta = real(y_z_theta_theta);
y_z_theta_theta(1) = y_z_theta_theta(1) + iv(-err,err);

out.y_phi_theta_theta = clip_tail(y_phi_theta_theta);
out.y_w_theta_theta = clip_tail(y_w_theta_theta);
out.y_u_theta_theta = clip_tail(y_u_theta_theta);
out.y_z_theta_theta = clip_tail(y_z_theta_theta);


