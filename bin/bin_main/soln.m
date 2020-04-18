function out=soln(x,s)
% out = soln(x,s)
%
% Returns the solution of bvp problem where the domain was split in half
%
% Input "x" is the value where the solution is evaluated and "s" is a
% stucture described in the STABLAB documenation

if x < 0
    x = s.side*s.I*(x/s.L);
    temp = deval(s.sol,x);
    out = temp(s.larray,:);
else
   x = s.side*s.I*(x/s.R);
   temp = deval(s.sol,x);
   out = temp(s.rarray,:);
end




