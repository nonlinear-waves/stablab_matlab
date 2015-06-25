function out = evan_root_bisector(fun,a,b,tol)
%
% out = evan_root_bisector(fun,a,b,tol)
%
% Returns the root of the Evans function with absolute tolerance tol
%
% Input is a left and right endpoint bounding the root, the tolernace, and
% a function handle that takes as input a and b and returns the values of
% the Evans function evaluated at the points a, b, and 0.5(a+b), as well
% as the points. [f,dom] = fun(a,b), where f = [D(a),D(0.5*(a+b)),D(b)],
% and dom = [a,0.5*(a+b),b];


% compute the Evans function on the points [a,0.5*(a+b),b].
[f,dom] = fun(a,b);
if f(1)*f(3) >= 0
    disp(dom);
    disp(f);
    error('Root is not bracketed.');
end

% update the interval based on compute values
if sign(f(1))==sign(f(2))
    a = dom(2);
elseif sign(f(2))==sign(f(3))
    b = dom(2);
else
    out = dom(2);
    return
end

% if width of interval is smaller than tolerance, return midpoint
if abs(a-b) < tol
   out = dom(2);
   return
end

% recursion step
out = evan_root_bisector(fun,a,b,tol);
