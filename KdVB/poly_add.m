function out = poly_add(a,b)

    % determine the size of the input and make the orientation standard
    [sx,sy] = size(a);
    if sy > sx
       a = a.'; 
    end
    if size(a,2) > 1
       error('vector input only'); 
    end
    
    [sx,sy] = size(b);
    if sy > sx
       b = b.'; 
    end
    if size(b,2) > 1
       error('vector input only'); 
    end
    
    

    N_a = length(a)-1; % degree of a poly
    N_b = length(b)-1; % degree of b poly

    
    % make the vectors the same size and add them
    if N_a < N_b
       a = [a;iv(zeros(N_b-N_a,1))];
    elseif N_b < N_a
       b = [b;iv(zeros(N_a-N_b,1))];  
    end
    
        
    out = a+b;