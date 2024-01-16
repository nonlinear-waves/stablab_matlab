function out = maxi(a,b)

% make sure that a,b are not NaN and take the maximum.

if sum(isnan(a))> 0
    disp(a)
   error('NaN encountered'); 
end

if nargin < 2  
    out = max(a);
else
    if sum(isnan(b))> 0
        disp(b)
        error('NaN encountered'); 
    end
    out = max(a,b);
end

