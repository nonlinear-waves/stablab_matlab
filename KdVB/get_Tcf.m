function Tcf = get_Tcf(M)

% Given an Mx1 vector, Tcf*v gives the coefficients of the
% polynomial interpolant of the function that takes on values v at the 
% Chebyshev nodes $x_j = (j+1/2)*pi/M$, $j = 0,...,M-1$.

% rigorous constants
half = iv(1)/2;
pie = iv('pi');

% Matrix for obtaining the Chebyshev coefficients
Id2 = (iv(2)/M)*speye(M);
Id2(1,1) = Id2(1,1)/2;
theta = ((0:1:M-1)+half)*pie/M;
Tcf = Id2*cos(theta.'*(0:1:M-1)).';