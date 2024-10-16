function [P,Q] = projection1(matrix,posneg,eps)
% [P,Q] = projection1(matrix,posneg,eps)
%
% Returns a projector P
%
% Input "matrix" is the matrix from which the eigenprojection comes,
% "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
% sought respectively. The input eps gives a bound on how small the eigenvalues sought
% can be, which is desirable when a zero mode should be avoided.

[R,D] = eig(matrix); 
L = inv(R);
P = zeros(size(R));

if posneg == 1
    index = find(real(diag(D))>eps).';
elseif posneg == -1
    index = find(real(diag(D))<eps).';
elseif posneg == 0
    index = find(abs(real(diag(D)))<eps).';
end

for j=index
    P = P + R(:,j)*L(j,:);
end

Q = P*R(:,index);

