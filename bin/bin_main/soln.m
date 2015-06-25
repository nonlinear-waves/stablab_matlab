function out=soln(x,s)
% out = soln(x,s)
%
% Returns the solution of bvp problem where the domain was split in half
%
% Input "x" is the value where the solution is evaluated and "s" is a
% stucture described in the STABLAB documenation

if x < 0
    x = s.side*s.I/s.L*x;
    temp = deval(s.sol,x);
    out = temp(s.larray,:);
else
   x = s.side*s.I/s.R*x;
   temp = deval(s.sol,x);
   out = temp(s.rarray,:);
end




