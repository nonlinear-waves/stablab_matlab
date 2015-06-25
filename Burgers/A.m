function out = A(x,lambda,s,p)
% function out = A(x,lambda,s,p)
%
% Evans matrix for Burgers system in unintegrated coordinates.

a=.5*(p.ul-p.ur);
cc=.5*(p.ul+p.ur); % wave speed
u=cc-a*tanh(a*x/2); % profile
uder=(-a^2/2)*sech(a*x/2)^2; % profile derivative

if strcmp(p.integrated,'on')
   out=[
            0                           1;
            lambda                 u-cc
       ];  
else
    out=[
                0                           1;
                lambda+uder        u-cc
           ];
end
   
