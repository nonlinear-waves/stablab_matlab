function [x_forward,x_backward] = get_conjugate_points(p,sol_mani,sol_mani_backward,muI,y0E0)

    % -------------------------------------------------------------------------
    % Compute $E_-^u$.
    % -------------------------------------------------------------------------

    options = odeset('RelTol',1e-12,'AbsTol',1e-18);

    %
    % backward part
    %

    % the ODE for the part backward from the first parametrization evaluation
    fun = @(x,y)(A_plucker(jac_OM(deval(sol_mani_backward,x),p))-muI)*y;

    % the solution to the ODE on the backward part
    sol_E0 = ode15s(fun,[sol_mani_backward.x(1),sol_mani_backward.x(end)],y0E0,options);

    % find indices that bound where a zero occurs
    ind = find(sol_E0.y(1,1:end-1).*sol_E0.y(1,2:end) < 0);
    
    % use the bisection method to find any zeros
    if isempty(ind)
        x_backward = [];
    else
        x_backward = zeros(1,length(ind));
        for j = 1:length(ind)
            xM = local_midpoint(sol_E0,sol_E0.x(ind(j)),sol_E0.x(ind(j)+1));
            x_backward(j) = xM;
        end  
    end

    %
    % forward part
    %

    % the ODE for the forward part from the first parametrization
    % evaluation
    fun = @(x,y)(A_plucker(jac_OM(deval(sol_mani,x),p))-muI)*y;
    
    % initialize the Plucker coordinates with the last values coming from
    % the backward computation
    y0E = deval(sol_E0,sol_E0.x(end));
    
    % compute the ODE slution for the forward part
    sol_E = ode15s(fun,[sol_mani.x(1),sol_mani.x(end)],y0E,options);

    % find the indices that bound the zeros
    ind = find(sol_E.y(1,1:end-1).*sol_E.y(1,2:end)< 0);

    % use the bisection method to find any zeros
    if isempty(ind)
        x_forward = [];
    else
        x_forward = zeros(1,length(ind));
        for j = 1:length(ind)
            xM = local_midpoint(sol_E,sol_E.x(ind(j)),sol_E.x(ind(j)+1));
            x_forward(j) = xM;
        end  
    end

end


%{
   find a zero between xL and xR using the bisection method
%}
function xM = local_midpoint(sol,xL,xR)
    
    % verify that there is a zero straddled
    fL = deval(sol,xL);
    fR = deval(sol,xR);
    if  fL(1)*fR(1) >= 0
       error('Zero is not straddled'); 
    end
    
    % use the bisection method to find the zero
    while xR-xL > 1e-11
       xM = (xL+xR)/2;
       temp = deval(sol,xM);
       fM = temp(1);
       if fR*fM > 0
           xR = xM;
       else
           xL = xM;
       end
    end

end
