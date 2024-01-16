function out = eval_cf_z(x,cf,a_x,b_x)

    % evaluate a Chebyshev polynomial corresponding to a function on [a_x,b_x]. 

    if size(x,2) > 1
       x = x.'; 
    end
    if size(cf,2) > 1
       cf = cf.'; 
    end
    N = length(cf);
    M=length(x);
    xtilde = (2*x-(a_x+b_x))/(b_x-a_x);
    for i=1:M
        xtilde(i)=intersect(xtilde(i),infsup(-1,1));
    end
    theta = acos(xtilde);
    
    T = cos(theta*(0:1:N-1));
    out = T*cf;
    
end
