
% Deval a function that is only defined on the left.  For example, if a 
% deval function is only defined from -10 to 0, and returns 8 values,
% this retruns the four values associated with any x.
function out = deval_left(eig, x)
    eigsize = size(deval(eig,0),1)/2;
    out = zeros(size(x,2),eigsize);
    L = 300;
    
    for i=1:size(x,2)
        
        % Right side
        if (x(1,i) > 0)
            out2 = deval(eig, x(1,i)-L);
            out(i,:) = [out2(5,1); ...
                        out2(7,1); ...
                        out2(6,1); ...
                        out2(8,1)];
            
        % Left side
        else 
            out2 = deval(eig, x(1,i));
            out(i,:) = [out2(1,1); ...
                        out2(3,1); ...
                        out2(2,1); ...
                        out2(4,1)];
        end
    end
end
