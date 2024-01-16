function [theta_R,r_R,sol_u] = fixed_x_R(nu,r_scale,x_R,theta_0,theta_R,r_R)

    % -------------------------------------------------------------------------
    % Find the appropriate values non-rigorously first
    % -------------------------------------------------------------------------

    % tolerance for solving the ODEs
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);

    fsolve_options.algorithm = optimset('Display','off', 'Jacobian','off', ...
            'Algorithm','Levenberg-Marquardt','TolFun',1e-14);

    bvp_options = bvpset('RelTol', 1e-13, 'AbsTol', 1e-13,'Nmax', 20000);
    mat = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
    
    % get the right manifold non-rigorously
    dr_nr = right_manifold_non_rigorous(mid(nu),r_scale);

    % get the left manifold non-rigorously
    dl_nr = left_manifold_single(mid(nu),theta_0,0);

    z0 = [theta_R,r_R];
  
    % solve for the profile from the left to the right non-rigorously
    sol_profile = ode15s(@(x,y)u_ode(x,y,mid(nu)),[0,x_R],...
        [dl_nr.phi_0;dl_nr.w_0;0;0],options);

    % solve for theta and r that correspond to a connection at the right
    [z,f,exitflag] = fsolve(@(u)fun_find_theta(u,dr_nr, ...
        sol_profile.y(1:2,end)),z0, ...
       fsolve_options.algorithm);
   
   

    % Choose theta \in [0,2\pi] and r>0
    if z(2) < 0
        z(2) = -z(2);
        z(1) = z(1)+pi;
    end

    while z(1) > 0
       z(1) = z(1)-2*pi; 
    end
    z(1) = z(1)+2*pi;

    yR = dr_nr(z(1),z(2));
    
    theta_R = z(1);
    r_R = z(2);
    

    %
    % solve the ricatti equation
    %

    data.phi0 = dl_nr.phi_0;
    data.w0 = dl_nr.w_0;
    data.u1 = yR(3);
    data.z1 = yR(4);
    bc = @(ya,yb)bc_u_ode(ya,yb,data);
    guess = @(x)mat*[deval(sol_profile,x);yR(3:4)];
    x_dom = linspace(0,x_R,30);
    solinit = bvpinit(x_dom,guess);
    sol_u = bvp5c(@(x,y)u_ode(x,y,mid(nu)),bc,solinit,bvp_options);



end




