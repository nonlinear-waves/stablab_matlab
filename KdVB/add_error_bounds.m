function poly_F = add_error_bounds(cf_F,delta_x,interpolation_error,truncation_error)
% Add the inerpolation error to the constant coefficient of the Chebyshev
% interpolant.
cf_F(1,:) = cf_F(1,:)+iv(-interpolation_error,interpolation_error);

%{
        Sum the truncated series, with error bound, to obtain an enclosure of
        the function F.
%}
poly_F = iv(zeros(size(cf_F,1),1));
for j = 1:fliplr(1:size(cf_F,2))
    poly_F = poly_F + delta_x^(j-1)*cf_F(:,j);
end

% Add the series truncation error.
poly_F(1) = poly_F(1)+iv(-truncation_error,truncation_error);
end