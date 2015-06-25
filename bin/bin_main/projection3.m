function [P,initial] = projection3(matrix,posneg,dud)
% [P,initializing vector] = projection3(matrix,posneg,eps)
%
% Returns a projector P for the single largest growth or decay mode.
% Used for compound Evans function method.
%
% Input "matrix" is the matrix from which the eigenprojection comes,
% "posneg" is 1,-1, or 0 if the unstable, stable, or center space is
% sought respectively.

[R,D] = eig(matrix); 
L = inv(R);

if posneg==1
    index = find(real(diag(D)) > max(real(diag(D)))-10^(-12));
elseif posneg==-1
    index = find(real(diag(D)) < min(real(diag(D)))+10^(-12));
end

P =  (R(:,index)*L(index,:))/(L(index,:)*R(:,index));
initial = R(:,index)/norm(R(:,index));

