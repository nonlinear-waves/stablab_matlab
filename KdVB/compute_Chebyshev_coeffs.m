function cf_F = compute_Chebyshev_coeffs(F,theta)
    [~,J]=size(F);
    %{
     Create the Chebyshev matrix operator Tcf that takes as input the
     evaluation of the function at the Chebyshev nodes and returns the
     coefficients of the Chebyshev interpolant.
    %}
    Id2 = (iv(2)/J)*speye(J);
    Id2(1,1) = Id2(1,1)/2;
    Tcf = Id2*cos(theta.'*(0:1:J-1)).';
    %    Compute the Chebyshev interpolant coefficients.
    cf_F = Tcf*F.';    
end