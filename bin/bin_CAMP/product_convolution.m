function out = product_convolution(varargin)
        % This code takes variable input and returns the convolution. For
        % example, product_convolution(u,v,w), where u, v, and w are m by
        % n matrices, returns 
        % $ \sum_{j=0}^m\sum_{k=0}^n \sum_{l=0}^j\sum_{s=0}^k
        % u_{m-j,n-k}v_{j-l,k-s}w_{l,s}$.

        m = length(varargin);
        if m == 2
            out = final_convolution(varargin{1},varargin{2});
        else
            S = intermediate_convolution(varargin{end-1},varargin{end});
            for j = 2:m-2
               S =  intermediate_convolution(varargin{end-j},S);
            end
            out = final_convolution(varargin{1},S);
        end
       
end

function sm = final_convolution(u,S)
   % takes as input u and S, which are m by n matrices, and returns
   % $ \sum_{j=0}^m\sum_{k=0}^n u_{m-j,n-k}S_{j,k}
   
    B = u.*rot90(S,2);
    sm = sum(sum(B));
    
end

function S = intermediate_convolution(v,w)
    % Takes as input v and w, which are m by n matrices, and returns an
    % m by n matrix whose (j,k) entry is $\sum_{l=0}^j\sum_{s=0}^k
    % v_{j-l,k-s}w_{l,s}.
    
    m = size(v,1)-1;
    n = size(v,2)-1;
    S = iv(zeros(m,n));
    for j = 0:m
        for k = 0:n
            S(j+1,k+1) = sum(sum(v(1:j+1,1:k+1).*rot90(w(1:j+1,1:k+1),2)));
        end
    end
    
end









