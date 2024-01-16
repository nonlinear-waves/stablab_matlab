function out = sin_of_poly(cf_poly)

%{
function out = cos_of_poly(cf_poly)

This function takes as input a vector of real-valued coefficients for 
a polynomial given in the Chebyshev basis. For example, if cf_poly = [1,2,3],
then the polynomial is 1*T0(x)+2*T1(x)+3*T2(x), where T0(x) = 1, T1(x) = x,
and T2(x) = 2x^2-1.

This function returns a Chebyshev polynomial enclosure of cos(cf_poly).

%}

% get a bound on the modulus of the interpolated value
bound = iv(sup(sum(abs(cf_poly))));

N = 25; % MUST BE ODD

% load the factorials
fac = iv(ones(N+1,1));
for k = 1:N
    fac(k+1) = fac(k)*k;
end

% get the truncation error bound for sine expansion: exp(bound)- (partial sum)
term = iv(0);
for k = fliplr(0:N)
    term = term+bound^k/fac(k+1);
end
truncation_error = exp(bound)-term;

% compute the powers of the polynomial input
x{N+1} = iv(1);
x{1} = iv(1);
for m = 1:N
    x{m+1} = poly_mult(x{m},cf_poly);
end

% partial sum of Taylor expansion of sine
out = iv(0);
for m = fliplr(0:1:floor(N/2)) % N must be even
    out = poly_add(out,(-1)^m*x{2*m+2}/fac(2*m+2));
end

% add the truncation error
out(1) = out(1)+truncation_error;

out = clip_tail(out);









