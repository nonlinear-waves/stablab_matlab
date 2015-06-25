function out = Akadj(x,lambda,s,p)
% Aadj(x,lambda,s,p)
%
% Returns the conjugate transpose of the matrix, s.A(x,lambda,s,p)

out = -s.Ak(x,lambda,s,p)';