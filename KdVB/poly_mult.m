function out = poly_mult(a,b)
% coefficients = poly_mult(a,b)
%
% Takes as input the coefficient vector a = [a_0, a_1,..., a_m] and the
% coefficient vector b = [b_0, b_1, ..., b_n] for Chebyshev polynomials
% a_0*T_0 + a_1*T_1 + ... + a_m*T_m and b_0*T_0 + b_1*T_1 + ... + b_n*T_n.
%
% Returns the coefficients for a Chebyshev polynomial enclosure of the
% product of the two Chebyshev polynomials. The domain of the
% polynomials is [-1,1]. 
% 
% For efficieny, the polynomial whose coefficients are returned is made of
% smaller degree with rigorous error of the truncation accounted for.

% ensure the coefficient vectors are column vectors.
if size(a,2) > size(a,1)
   a = a.'; 
end
if size(b,2) > size(b,1)
   b = b.'; 
end

% make sure the input is that of a one-dimensional polynomial
if size(a,2) > 1
   error('input is not valid'); 
end
if size(b,2) > 1
   error('input is not valid'); 
end

% Obtain the degree of the polynomails given.
deg_a = length(a)-1;
deg_b = length(b)-1;

% make the polynomials have the same degree
if deg_a < deg_b
   a = [a;iv(zeros(deg_b-deg_a,1))];
elseif deg_b < deg_a
   b = [b;iv(zeros(deg_a-deg_b,1))];  
end

% degree of the polynomials
deg = maxi(deg_a,deg_b);

% compute the outer product of a and b, divided by 2
outer = (a/2)*b.';

% sum the terms for $T_{i+j}$ 
out = sum(spdiags(rot90(outer),-deg:deg));

% sum the terms for $T_{|i-j|} = 0$.
diag_sum = sum(spdiags(outer,-deg:deg));

% update the coefficient for $T_0$.
out(1) = out(1)+diag_sum(deg+1);

% update the coefficients for $T_1$ up to $T_{deg}$.
out(2:deg+1) = out(2:deg+1)+ fliplr(diag_sum(1:deg)+fliplr(diag_sum(deg+2:end)));

% return a column vector.
out = out.';

% Find out how many terms of the tail of the polynomial are small enough
% to clip with an added error.
out = clip_tail(out);




