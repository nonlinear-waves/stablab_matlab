function out = clip_tail(out)

if length(out) > 2
    
    % find the maximum value of the input vector
    mx = maxi(sup(abs(out)));
    remove_index = 0;
    err = 0;
    while err <= 1e-16*mx
        err = err + sup(abs(out(end-remove_index)));
        remove_index = remove_index+1;
        if remove_index == length(out)-1
           break 
        end
    end
    remove_index = maxi(0,remove_index-1);

    % figure out the truncation error
    truncation_error = sum(abs(out(end-remove_index+1:end)));


    % add the truncation error to the polynomial.
    if maxi(sup(abs(imag(out)))) > 0
        out(1) = out(1)+(1+1i)*iv(-truncation_error,truncation_error);
    else
        out(1) = out(1)+iv(-truncation_error,truncation_error);
    end
    
    % remove the now extra terms
    out = out(1:end-remove_index); 
   


end