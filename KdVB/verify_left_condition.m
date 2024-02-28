function verify_left_condition(theta_0,mu_left,mu_right,phi_poly,phi_x_poly, ...
    u_poly,u_x_poly,r_newton,nu_left,nu_right)



        % add on error to the profiles evalauted at x = 0
        phi_poly(1) = phi_poly(1) + iv(-r_newton,r_newton);
        phi_x_poly(1) = phi_x_poly(1) + iv(-r_newton,r_newton);
        u_poly(1) = u_poly(1) + iv(-r_newton,r_newton);
        u_x_poly(1) = u_x_poly(1) + iv(-r_newton,r_newton);

        % get grid of intervals for evaluating the profiles
        theta_interp = linspace(0,sup(2*iv('pi')),1001);
        theta_interp = iv(theta_interp(1:end-1),theta_interp(2:end)).';

        % evaluate the profiles at the grid points
        phi = cos(theta_interp*(0:1:length(phi_poly)-1))*phi_poly;
        phi_x = cos(theta_interp*(0:1:length(phi_x_poly)-1))*phi_x_poly;
        u = cos(theta_interp*(0:1:length(u_poly)-1))*u_poly;
        u_x = cos(theta_interp*(0:1:length(u_x_poly)-1))*u_x_poly;

        % get the range of nu values
        nu_range = iv(nu_left,nu_right);

        % condition > 0 indicates success
        condition = 1;

        %
        % check the first condition
        %

        mu_range = iv(mu_left,mu_right);

        str = '';
        if sup(abs(theta_0/(1-theta_0/(6-4*mu_range)))) >= 2

            str = 'condition 1';
            condition = -1;
        end

        %
        % check the second condition
        %

        phi_xx_times_nu = -phi_x+(phi.^2-1)/2;


        if max(sup(phi_xx_times_nu)) >= 0
            str = [str,' condition 2'];
            condition = -1;
        end

        %
        % check the third condition
        %

        v = u_x./u;

        if sum(isnan(mid(v))) > 0
            error('failed: divide by zero?');
        end

        max_sup = max(sup(v));
        
        if max_sup >= 0
            str = [str, 'condition 3'];
            condition = -1;
        end


        %
        % check the fourth condition
        %

        term = -atan(-sqrt(2./(-phi_x)).*v)+ ...
            (1/sqrt(mu_range))*log((1+sqrt(theta_0/2))/(1-sqrt(theta_0/2)));

        if sum(isnan(mid(term))) > 0
            error('failed: divide by zero?');
        end

        max_term = max(sup(term));

        if max_term >= 0

            str = [str, ' condition 4'];
            condition = -1;
        end

        %
        % check if at least one of the conditions failed
        %

        if condition == -1
            disp(str);
            error('Failed to verify no conjugate points to the left of x = 0.')

        end




