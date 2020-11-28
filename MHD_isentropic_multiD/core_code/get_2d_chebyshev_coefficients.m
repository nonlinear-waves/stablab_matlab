function [cf,fun] = get_2d_chebyshev_coefficients(f,a_x,b_x,a_y,b_y,opt)

sz = size(f);
if ~(length(sz)==2)
   error('input must be a matrix'); 
end

if strcmp(opt,'first kind')
    
    N_x = size(f,1);
    N_y = size(f,2);
    
    theta_x = pi*(1:2:2*N_x-1)/(2*N_x);
    theta_y = pi*(1:2:2*N_y-1)/(2*N_y);


    
    Tx = cos((0:1:N_x-1).'*theta_x)
    Ty = cos(theta_y.'*(0:1:N_y-1))

    cf = (4/(N_x*N_y))*Tx*f*Ty;
    cf(1,:) = cf(1,:)/2;
    cf(:,1) = cf(:,1)/2;

    

    fun = @(x,y)eval_cf(x,y,cf,a_x,b_x,a_y,b_y);
    
elseif strcmp(opt,'FirstKindExtrema')
    
    
    
end


% Transformation to get Chebyshev coefficients

function out = eval_cf(x,y,cf,a_x,b_x,a_y,b_y)

    
    N_x = size(cf,1);
    N_y = size(cf,2);

    xtilde = (x-0.5*(a_x+b_x))/(0.5*(b_x-a_x));
    ytilde = (y-0.5*(a_y+b_y))/(0.5*(b_y-a_y));
    
    theta_x = acos(xtilde);
    theta_y = acos(ytilde);
    
    Tx = cos(theta_x.'*(0:1:N_x-1));
    Ty = cos((0:1:N_y-1).'*theta_y);
    
    
    out = Tx*cf*Ty;

    







