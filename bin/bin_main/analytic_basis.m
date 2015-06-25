function [out, projects] = analytic_basis(projection,x,preimage,s,p,A,posneg,eps,Q,p_old_in)
% [out, projects] = analytic_basis(projection,x,preimage,s,p,A,posneg,eps,Q,p_old_in)
%
% Returns an analytic basis at specified infinity using the method of Kato.
%
% Input "projection" is a function handle to the projection function to be
% used, "x" is the numerical value of infinity, "preimage" is the contour
% on which the Evans function is computed, "s" and "p" are structures
% explained in the STABLAB documentation, "A" is a function handle to the
% Evans matrix, "posneg" is 1 or -1 determining which space the
% projection function should return, and "eps" is the tolerance in the
% projection function. If input Q and p_old_in are specified, then the
% first value of the analtyic basis is set to Q with projection p_old_in.
% This allows for continuing or filling in a previous analtyic basis
% computation. When these two inputs are not specified, they are chosen
% automtically.

iterations = size(preimage,2);
[p_old, Q1] = projection(A(x,preimage(1),s,p),posneg,eps);
[n,k]=size(Q1);
out = zeros(n,k,iterations);
projects = zeros(size(p_old,1),size(p_old,2),iterations);


if nargin < 9
    out(:,:,1) = Q1;
    projects(:,:,1) = p_old;
else
    out(:,:,1) = Q;
    projects(:,:,1) = p_old_in;
    p_old = p_old_in;
end

for j=2:iterations
    proj = projection(A(x,preimage(j),s,p),posneg,eps);
%     out(:,:,j) = proj * (eye(n) + proj * p_old - p_old * proj) * out(:,:,j-1);
    out(:,:,j) = proj * (eye(n) + 0.5 * p_old * (eye(n)- proj)) * out(:,:,j-1);
    projects(:,:,j) = proj;
    p_old = proj;
end
 




