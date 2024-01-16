function out = iv(num1,num2)
% function out = iv(num)
%
% returns a mathematically rigorous enclosure of the number of interval
% given.

if nargin == 1
    out = intval(num1);
elseif nargin == 2
    out = hull(intval(num1),intval(num2));
end