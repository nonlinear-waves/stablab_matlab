function [P,Q1] = projection2(matrix,posneg,eps)
% [P,Q1] = projection2(matrix,posneg,eps)
%
% Returns a projector P and spanning set Q1 of the invariant subspace
% associated with the given matrix and specified subspace.
%
% Input "matrix" is the matrix from which the eigenprojection comes,
% "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
% sought. The input eps gives a bound on how small the eigenvalues sought
% can be, which is desirable when a zero mode should be avoided.

% Uses Schur decomposition to get a basis for the generalized eigenspace

[U,T] = schur(matrix,'complex');
E = ordeig(T);
k = length(find(posneg*real(E)>eps));
US = ordschur(U,T,posneg*real(E)>eps);
Q1 = US(:,1:k);

[U,T] = schur(-matrix,'complex');
E = ordeig(T);
k = length(find(posneg*real(E)>-eps));
US = ordschur(U,T,posneg*real(E)>-eps);
Q2 = US(:,1:k);

R = [Q1 Q2];
L = inv(R);

P = zeros(size(matrix));
for k=1:size(Q1,2)
    P = P + R(:,k)*L(k,:);
end


