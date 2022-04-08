function [P,Q] = projection5(matrix,posneg,eps)
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

e.kl=2;
e.kr=1;

kl=e.kl;
kr=e.kr;
 
eg = diag(D);
[unused,ind] = sort(real(eg));
 
if posneg == 1
   % left side
   ind = flipud(ind);
   index = ind(1:kl);
%    index = ind([2,3]);
elseif posneg == -1
   % positive side
%    index = ind(1);
   index = ind(1:kr);
end
 
for j=index
    P = P + R(:,j)*L(j,:);
end
 
Q = P*R(:,index);
