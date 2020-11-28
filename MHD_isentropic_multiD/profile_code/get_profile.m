function [s,p] = get_profile(p,end_point_tol)

% profile
s.something = 0;
num_inf = 2;
if nargin < 2
    end_point_tol = 1e-8;
end
max_err = 1+end_point_tol;
while max_err > end_point_tol

    [s,p] = profile_solve_pseudo(s,p,num_inf);

    max_err = abs(s.sol.deval(s.L)-1);
    num_inf = 2*num_inf;
    if num_inf > 128
       error('failed to solve the profile'); 
    end
end