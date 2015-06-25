function out = Aadj(x,lambda,s,p)
% Aadj(x,lambda,s,p)
%
% Returns the conjugate transpose of the matrix, s.A(x,lambda,s,p)

out = -s.A(x,lambda,s,p)';